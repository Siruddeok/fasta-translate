"""FASTA 서열을 검증하고 아미노산으로 번역해 CSV로 저장한다.

사용 예:
    python fasta_translate.py input.fasta -o result.csv
    python fasta_translate.py mrna.fasta --start first-atg
    python fasta_translate.py genome_frag.fasta --start longest-orf --table 11 --alt-starts
    python fasta_translate.py input.fasta --strict --table 2

번역 시작 위치 (--start)
    frame1      : 1번 염기부터 서열 끝까지 번역 (CDS를 넣는 경우, 기본값)
    first-atg   : 첫 개시 코돈부터 그 틀의 첫 종결 코돈까지 번역 (5'UTR이 있는 mRNA)
    longest-orf : 정방향 3개 틀에서 가장 긴 ORF를 찾아 번역
    first-atg/longest-orf는 기본적으로 ATG만 개시 코돈으로 찾는다.
    --alt-starts를 주면 --table 유전 코드가 허용하는 모든 개시 코돈(GTG, TTG 등)을 찾는다.
"""

import argparse
import csv
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError, unambiguous_dna_by_id
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

DNA_BASES = set("ATGC")
AMBIGUOUS = set("N")
GAP_CHARS = set("-.")
START_MODES = ("frame1", "first-atg", "longest-orf")

COLUMNS = [
    "ID",
    "Description",
    "Seq_Type",
    "Length",
    "GC_Content",
    "Start_Pos",
    "End_Pos",
    "Is_Valid",
    "Errors",
    "Warnings",
    "Protein_Full",
    "Protein_To_Stop",
    "Protein_Length",
]


class GeneticCode:
    """--table 번호에 맞는 개시/종결 코돈 목록."""

    def __init__(self, table_id):
        if table_id not in unambiguous_dna_by_id:
            raise ValueError(f"지원하지 않는 유전 코드 번호: {table_id}")
        table = unambiguous_dna_by_id[table_id]
        self.id = table_id
        self.starts = set(table.start_codons)
        self.stops = set(table.stop_codons)
        self.starts_label = "/".join(
            ["ATG"] + sorted(c for c in self.starts if c != "ATG")
        )


def normalize(seq):
    """갭 제거, DNA/RNA 판별, 허용 문자 검사.

    반환: (DNA로 정규화된 서열, 서열 종류, 오류 목록, 경고 목록)
    """
    errors = []
    warnings = []

    if any(c in GAP_CHARS for c in seq):
        seq = "".join(c for c in seq if c not in GAP_CHARS)
        warnings.append("갭 문자 제거됨 (좌표는 갭 제거 후 기준)")

    if not seq:
        return seq, "UNKNOWN", ["빈 서열"], warnings

    has_t = "T" in seq
    has_u = "U" in seq
    if has_t and has_u:
        errors.append("T와 U가 함께 존재")
        seq_type = "MIXED"
    elif has_u:
        seq_type = "RNA"
        seq = seq.replace("U", "T")
    else:
        seq_type = "DNA"

    allowed = DNA_BASES | AMBIGUOUS | ({"U"} if seq_type == "MIXED" else set())
    invalid = set(seq) - allowed
    if invalid:
        errors.append(f"허용되지 않은 문자: {''.join(sorted(invalid))}")

    n_count = seq.count("N")
    if n_count:
        warnings.append(f"모호 염기 N {n_count}개 (X로 번역됨)")

    return seq, seq_type, errors, warnings


def first_stop(seq, start, code):
    """start 위치 틀에서 첫 종결 코돈의 시작 인덱스. 없으면 None."""
    for i in range(start, len(seq) - 2, 3):
        if seq[i:i + 3] in code.stops:
            return i
    return None


def find_orf(seq, mode, code, alt_starts):
    """번역 구간 (start, end)을 0-based 반열린 구간으로 반환. 못 찾으면 None.

    frame1이면 서열 전체, 그 외에는 개시 코돈부터 첫 종결 코돈(포함)까지.
    종결 코돈이 없으면 서열 끝까지.
    """
    if mode == "frame1":
        return 0, len(seq)

    starts = code.starts if alt_starts else {"ATG"}

    def orf_end(start):
        stop = first_stop(seq, start, code)
        return len(seq) if stop is None else stop + 3

    if mode == "first-atg":
        for i in range(len(seq) - 2):
            if seq[i:i + 3] in starts:
                return i, orf_end(i)
        return None

    # longest-orf: 정방향 3개 틀, 각 ORF의 가장 앞 개시 코돈 기준
    best = None
    for frame in range(3):
        i = frame
        while i <= len(seq) - 3:
            if seq[i:i + 3] in starts:
                end = orf_end(i)
                codons = (end - i) // 3
                if best is None or codons > (best[1] - best[0]) // 3:
                    best = (i, end)
                i = end  # 같은 ORF 안의 뒤쪽 ATG는 더 짧으므로 건너뜀
            else:
                i += 3
    return best


def check_structure(seq, start, end, mode, code):
    """번역 구간의 구조를 검사한다. 반환: (오류 목록, 경고 목록)"""
    errors = []
    warnings = []
    coding = seq[start:end]
    remainder = len(coding) % 3
    trimmed_len = len(coding) - remainder

    if trimmed_len < 3:
        return ["번역할 코돈이 없음 (3염기 미만)"], warnings

    if remainder:
        warnings.append(f"길이가 3의 배수가 아님 (끝 {remainder}염기 버림)")

    first = coding[:3]
    if first not in code.starts:
        warnings.append(f"개시 코돈({code.starts_label}) 없음")
    elif first != "ATG":
        warnings.append(f"대체 개시 코돈 {first}로 시작")

    if first in code.stops:
        errors.append("첫 코돈이 종결 코돈 (번역 결과 없음)")
        return errors, warnings

    last = coding[trimmed_len - 3:trimmed_len]
    if last not in code.stops:
        if mode == "frame1":
            warnings.append("종결 코돈으로 끝나지 않음")
        else:
            warnings.append("ORF 안에 종결 코돈 없음 (서열 끝까지 번역)")

    if mode == "frame1":
        internal = [
            str(start + i + 1)
            for i in range(0, trimmed_len - 3, 3)
            if coding[i:i + 3] in code.stops
        ]
        if internal:
            warnings.append(f"중간 종결 코돈 위치(nt): {', '.join(internal)}")

    return errors, warnings


def translate(coding, code, start_is_chosen):
    """번역 구간을 번역한다. (전체 번역, 첫 종결 코돈 전까지 번역)을 반환.

    첫 코돈이 유전 코드가 허용하는 개시 코돈이면 Met로 번역한다.
      - start_is_chosen=True (first-atg/longest-orf): 개시 코돈으로 골랐으므로 항상
      - frame1: 완전한 CDS(종결 코돈 끝, 3의 배수, 중간 종결 없음)일 때만.
        CDS 중간 조각이 우연히 CTG 등으로 시작해도 Met로 바꾸지 않기 위함.
    """
    coding = coding[: len(coding) - len(coding) % 3]
    seq_obj = Seq(coding)
    try:
        protein = str(seq_obj.translate(table=code.id, cds=True))
        return protein + "*", protein
    except TranslationError:
        pass
    full = str(seq_obj.translate(table=code.id))
    to_stop = str(seq_obj.translate(table=code.id, to_stop=True))
    if start_is_chosen and coding[:3] in code.starts:
        full = "M" + full[1:]
        to_stop = "M" + to_stop[1:]
    return full, to_stop


def process_record(record, code, mode, alt_starts, strict, seen_ids):
    raw = str(record.seq).upper()
    seq, seq_type, errors, warnings = normalize(raw)

    if record.id in seen_ids:
        warnings.append("중복 ID")
    seen_ids.add(record.id)

    start = end = None
    if not errors:
        orf = find_orf(seq, mode, code, alt_starts)
        if orf is None:
            label = code.starts_label if alt_starts else "ATG"
            errors.append(f"개시 코돈({label})을 찾지 못함")
        else:
            start, end = orf
            struct_errors, struct_warnings = check_structure(seq, start, end, mode, code)
            errors += struct_errors
            warnings += struct_warnings

    if strict and warnings:
        errors.extend(warnings)
        warnings = []

    if errors:
        full, to_stop = "", ""
    else:
        full, to_stop = translate(seq[start:end], code, start_is_chosen=mode != "frame1")

    return {
        "ID": record.id,
        "Description": record.description,
        "Seq_Type": seq_type,
        "Length": len(raw),
        "GC_Content": round(gc_fraction(seq), 4) if seq else "",
        "Start_Pos": "" if start is None else start + 1,
        "End_Pos": "" if end is None else end,
        "Is_Valid": not errors,
        "Errors": "; ".join(errors),
        "Warnings": "; ".join(warnings),
        "Protein_Full": full,
        "Protein_To_Stop": to_stop,
        "Protein_Length": len(to_stop),
    }


def write_csv(rows, path):
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def parse_args():
    parser = argparse.ArgumentParser(
        description="FASTA 서열 검증 후 아미노산 번역 결과를 CSV로 저장"
    )
    parser.add_argument("input", help="입력 FASTA 파일 경로")
    parser.add_argument(
        "-o", "--output", help="출력 CSV 경로 (기본: <입력파일명>_translated.csv)"
    )
    parser.add_argument(
        "--table", type=int, default=1, help="NCBI 유전 코드 번호 (기본 1: 표준)"
    )
    parser.add_argument(
        "--start", choices=START_MODES, default="frame1",
        help="번역 시작 위치 (기본 frame1: 1번 염기부터)",
    )
    parser.add_argument(
        "--alt-starts", action="store_true",
        help="first-atg/longest-orf에서 유전 코드의 대체 개시 코돈도 찾기",
    )
    parser.add_argument(
        "--strict", action="store_true", help="경고도 오류로 취급 (완전한 CDS만 유효)"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input)
    if not input_path.is_file():
        sys.exit(f"오류: 파일을 찾을 수 없음: {input_path}")

    try:
        code = GeneticCode(args.table)
    except ValueError as e:
        sys.exit(f"오류: {e}")

    output_path = (
        Path(args.output)
        if args.output
        else input_path.with_name(f"{input_path.stem}_translated.csv")
    )

    rows = []
    seen_ids = set()
    try:
        for record in SeqIO.parse(input_path, "fasta"):
            rows.append(
                process_record(record, code, args.start, args.alt_starts, args.strict, seen_ids)
            )
    except (ValueError, UnicodeDecodeError) as e:
        sys.exit(f"오류: FASTA 파일을 읽을 수 없음: {e}")

    if not rows:
        sys.exit("오류: FASTA 레코드가 없음 (빈 파일 또는 형식 오류)")

    write_csv(rows, output_path)

    valid = sum(r["Is_Valid"] for r in rows)
    with_warn = sum(1 for r in rows if r["Warnings"])
    print(f"전체 {len(rows)}개 / 유효 {valid}개 / 무효 {len(rows) - valid}개 / 경고 포함 {with_warn}개")
    print(f"결과 저장: {output_path}")


if __name__ == "__main__":
    main()
