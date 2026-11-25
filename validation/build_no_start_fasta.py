"""개시 코돈(ATG)으로 시작하지 않는 서열 검증용 FASTA를 만든다.

build_indel_fasta.py가 내려받은 RefSeq 캐시(data/*.gb)를 사용한다.
사용: python validation/build_no_start_fasta.py
"""

from pathlib import Path

from Bio import SeqIO

HERE = Path(__file__).resolve().parent


def load(acc):
    record = SeqIO.read(HERE / "data" / f"{acc}.gb", "genbank")
    cds = next(f for f in record.features if f.type == "CDS")
    return record, cds


def main():
    hbb, hbb_cds = load("NM_000518")
    cftr, cftr_cds = load("NM_000492")
    mrna = str(hbb.seq).upper()
    cds = str(hbb_cds.extract(hbb.seq)).upper()
    cds_start = int(hbb_cds.location.start)          # 0-based, 5'UTR 길이
    cftr_seq = str(cftr_cds.extract(cftr.seq)).upper()

    cases = [
        ("A_HBB_mRNA_with_5UTR",
         f"HBB 전체 mRNA, 5'UTR {cds_start} nt 포함 (CDS는 {cds_start + 1}번부터)", mrna),
        ("B_HBB_CDS_no_ATG",
         "HBB CDS에서 ATG만 제거 (틀 유지)", cds[3:]),
        ("C_HBB_CDS_minus1nt",
         "HBB CDS 앞 1 nt 제거 (틀 어긋남)", cds[1:]),
        ("D_CFTR_internal_fragment",
         "CFTR CDS c.1501-1650 내부 조각 (코돈 경계에서 시작)", cftr_seq[1500:1650]),
        ("E_HBB_CDS_lowercase",
         "HBB CDS 소문자 (atg 시작)", cds.lower()),
        ("F_HBB_CDS_leading_N",
         "HBB CDS 앞에 NNN", "NNN" + cds),
        ("G_stop_first",
         "첫 코돈이 종결 코돈 + HBB CDS", "TAA" + cds),
        ("H_HBB_CDS_GTG_start",
         "HBB CDS의 ATG를 GTG로 (세균 대체 개시 코돈 형태)", "GTG" + cds[3:]),
    ]

    out = HERE / "no_start_test.fasta"
    with open(out, "w", encoding="utf-8") as f:
        for seq_id, desc, seq in cases:
            f.write(f">{seq_id} {desc}\n{seq}\n")
    print(f"HBB 5'UTR {cds_start} nt, CDS {len(cds)} nt / {len(cases)}개 서열 → {out}")


if __name__ == "__main__":
    main()
