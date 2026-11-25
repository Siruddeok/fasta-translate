"""공개된 실제 indel 돌연변이로 검증용 FASTA를 만든다.

1. NCBI RefSeq mRNA(GenBank)를 내려받아 CDS를 추출한다.
2. 문헌/ClinVar에 보고된 HGVS c. 표기대로 indel을 적용한다.
   적용 전에 RefSeq의 해당 위치 염기가 HGVS에 적힌 염기와 같은지 확인한다.
3. 야생형 + 돌연변이 서열을 indel_test.fasta로,
   HGVS p. 표기에서 계산한 기대 단백질 길이를 expected.csv로 저장한다.

사용: python validation/build_indel_fasta.py
"""

import csv
import io
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

from Bio import SeqIO

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# start/end: HGVS c. 좌표 (A of ATG = 1)
# ref: 삭제·중복되는 염기 (RefSeq와 대조). None이면 대조 없이 길이만 사용
# exp_wt / exp_mut: 종결 코돈 제외 단백질 길이(aa, Met 포함)
#   - frameshift p.XxxNYyyfs*K → 새 종결 코돈은 N번째 코돈에서 K번째(첫 변경 aa = 1) → 종결 직전까지 길이 N+K-2
MUTATIONS = [
    {
        "gene": "HBB", "acc": "NM_000518", "kind": "del", "start": 126, "end": 129, "ref": "CTTT",
        "hgvs_c": "c.126_129del", "hgvs_p": "p.Phe42Leufs*19",
        "note": "베타지중해빈혈, codon 41/42 (-TTCT)", "exp_wt": 147, "exp_mut": 42 + 19 - 2,
    },
    {
        "gene": "CFTR", "acc": "NM_000492", "kind": "del", "start": 1521, "end": 1523, "ref": "CTT",
        "hgvs_c": "c.1521_1523del", "hgvs_p": "p.Phe508del",
        "note": "낭포성 섬유증, ΔF508 (in-frame)", "exp_wt": 1480, "exp_mut": 1479,
    },
    {
        "gene": "BRCA1", "acc": "NM_007294", "kind": "del", "start": 68, "end": 69, "ref": "AG",
        "hgvs_c": "c.68_69del", "hgvs_p": "p.Glu23Valfs*17",
        "note": "유전성 유방암/난소암, 185delAG", "exp_wt": 1863, "exp_mut": 23 + 17 - 2,
    },
    {
        "gene": "BRCA2", "acc": "NM_000059", "kind": "del", "start": 5946, "end": 5946, "ref": "T",
        "hgvs_c": "c.5946del", "hgvs_p": "p.Ser1982Argfs*22",
        "note": "유전성 유방암/난소암, 6174delT", "exp_wt": 3418, "exp_mut": 1982 + 22 - 2,
    },
    {
        "gene": "HEXA", "acc": "NM_000520", "kind": "dup", "start": 1274, "end": 1277, "ref": "TATC",
        "hgvs_c": "c.1274_1277dup", "hgvs_p": "p.Tyr427Ilefs*5",
        "note": "테이-삭스병, 1278insTATC", "exp_wt": 529, "exp_mut": 427 + 5 - 2,
    },
    {
        "gene": "CCR5", "acc": "NM_000579", "kind": "del", "start": 554, "end": 585, "ref": None,
        "hgvs_c": "c.554_585del", "hgvs_p": "p.Ser185Ilefs*32",
        "note": "HIV 저항성, CCR5-Δ32 (32bp 결실)", "exp_wt": 352, "exp_mut": 185 + 32 - 2,
    },
]


def fetch_genbank(acc):
    """RefSeq GenBank 레코드를 내려받는다 (data/ 에 캐시)."""
    DATA_DIR.mkdir(exist_ok=True)
    cache = DATA_DIR / f"{acc}.gb"
    if cache.exists():
        return SeqIO.read(cache, "genbank")

    params = urllib.parse.urlencode(
        {"db": "nuccore", "id": acc, "rettype": "gb", "retmode": "text",
         "tool": "fasta_translate_validation"}
    )
    with urllib.request.urlopen(f"{EFETCH}?{params}", timeout=60) as resp:
        text = resp.read().decode("utf-8")
    cache.write_text(text, encoding="utf-8")
    time.sleep(0.4)  # NCBI 권장: API 키 없이 초당 3회 이하
    return SeqIO.read(io.StringIO(text), "genbank")


def extract_cds(record):
    cds_features = [f for f in record.features if f.type == "CDS"]
    if len(cds_features) != 1:
        raise ValueError(f"{record.id}: CDS feature가 {len(cds_features)}개")
    return str(cds_features[0].extract(record.seq)).upper()


def apply_indel(cds, m):
    s, e = m["start"], m["end"]
    target = cds[s - 1:e]
    if m["ref"] is not None and target != m["ref"]:
        raise ValueError(
            f"{m['gene']} {m['hgvs_c']}: RefSeq 염기 {target} ≠ HGVS 염기 {m['ref']}"
        )
    if m["kind"] == "del":
        return cds[:s - 1] + cds[e:], target
    if m["kind"] == "dup":
        return cds[:e] + target + cds[e:], target
    raise ValueError(f"알 수 없는 변이 종류: {m['kind']}")


def main():
    fasta_lines = []
    expected = []
    for m in MUTATIONS:
        record = fetch_genbank(m["acc"])
        cds = extract_cds(record)
        mutant, target = apply_indel(cds, m)
        diff = len(mutant) - len(cds)
        frame = "in-frame" if diff % 3 == 0 else "frameshift"
        print(f"{m['gene']:6} {record.id:13} CDS {len(cds):5} nt | {m['hgvs_c']:15} "
              f"{m['kind']} {target} ({diff:+d} nt, {frame})")

        wt_id, mut_id = f"{m['gene']}_WT", f"{m['gene']}_{m['hgvs_c']}"
        fasta_lines += [f">{wt_id} {record.id} wild-type CDS", cds]
        fasta_lines += [f">{mut_id} {record.id}:{m['hgvs_c']} {m['hgvs_p']} {m['note']}", mutant]
        expected.append({"ID": wt_id, "Accession": record.id, "HGVS_c": "", "HGVS_p": "",
                         "Frame": "", "Expected_Protein_Length": m["exp_wt"]})
        expected.append({"ID": mut_id, "Accession": record.id, "HGVS_c": m["hgvs_c"],
                         "HGVS_p": m["hgvs_p"], "Frame": frame,
                         "Expected_Protein_Length": m["exp_mut"]})

    (HERE / "indel_test.fasta").write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")
    with open(HERE / "expected.csv", "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(expected[0]))
        writer.writeheader()
        writer.writerows(expected)
    print(f"\n{len(expected)}개 서열 → {HERE / 'indel_test.fasta'}")


if __name__ == "__main__":
    try:
        main()
    except ValueError as e:
        sys.exit(f"오류: {e}")
