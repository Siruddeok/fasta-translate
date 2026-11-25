"""fasta_translate.py 결과를 HGVS 기반 기대값과 비교한다.

검사 항목
  1. 야생형: 유효, 경고 없음, 단백질 길이 = RefSeq 단백질 길이
  2. 돌연변이: 유효(번역 수행), 단백질 길이 = HGVS p. 표기에서 계산한 길이
  3. frameshift: 길이가 3의 배수가 아님/중간 종결 코돈 경고가 있어야 함
  4. 야생형 대비 단백질 서열: 변이 지점 앞부분은 동일해야 함
     (in-frame 결실은 해당 아미노산만 빠져야 함)

사용: python validation/compare_results.py
"""

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(path):
    with open(path, encoding="utf-8-sig", newline="") as f:
        return {row["ID"]: row for row in csv.DictReader(f)}


def check(row, exp, wt_row):
    problems = []
    length = int(row["Protein_Length"])
    if row["Is_Valid"] != "True":
        problems.append(f"무효 처리됨: {row['Errors']}")
    if length != int(exp["Expected_Protein_Length"]):
        problems.append(f"단백질 길이 {length} ≠ 기대 {exp['Expected_Protein_Length']}")

    if wt_row is None:  # 야생형
        if row["Warnings"]:
            problems.append(f"야생형인데 경고 있음: {row['Warnings']}")
        return problems

    # 돌연변이: 변이 지점(p. 표기의 첫 번호) 앞까지 야생형과 동일해야 함
    pos = int(re.search(r"\d+", exp["HGVS_p"]).group())
    wt_prot, mut_prot = wt_row["Protein_To_Stop"], row["Protein_To_Stop"]
    if mut_prot[:pos - 1] != wt_prot[:pos - 1]:
        problems.append(f"{pos}번 이전 아미노산이 야생형과 다름")

    if exp["Frame"] == "frameshift":
        if "3의 배수가 아님" not in row["Warnings"] and "중간 종결 코돈" not in row["Warnings"]:
            problems.append("frameshift인데 관련 경고 없음")
        new_aa = mut_prot[pos - 1] if len(mut_prot) >= pos else "-"
        expected_aa = exp["HGVS_p"][len(f"p.Xxx{pos}"):][:3]
        three_to_one = {"Leu": "L", "Val": "V", "Arg": "R", "Ile": "I"}
        if three_to_one.get(expected_aa) != new_aa:
            problems.append(f"{pos}번 아미노산 {new_aa} ≠ HGVS {expected_aa}")
    else:  # in-frame 결실: 해당 아미노산 하나만 빠짐
        if wt_prot[:pos - 1] + wt_prot[pos:] != mut_prot:
            problems.append("in-frame 결실 외 다른 차이 존재")
    return problems


def main():
    results = load(HERE / "indel_translated.csv")
    expected = load(HERE / "expected.csv")

    failures = 0
    print(f"{'ID':24} {'유효':4} {'길이':>5} {'기대':>5}  판정  경고")
    for seq_id, exp in expected.items():
        row = results.get(seq_id)
        if row is None:
            print(f"{seq_id:24} 결과 없음")
            failures += 1
            continue
        gene = seq_id.split("_")[0]
        wt_row = None if seq_id.endswith("_WT") else results.get(f"{gene}_WT")
        problems = check(row, exp, wt_row)
        failures += bool(problems)
        verdict = "PASS" if not problems else "FAIL"
        warn = row["Warnings"] or "-"
        print(f"{seq_id:24} {row['Is_Valid']:4} {row['Protein_Length']:>5} "
              f"{exp['Expected_Protein_Length']:>5}  {verdict}  {warn[:70]}")
        for p in problems:
            print(f"    └ {p}")

    print(f"\n{len(expected) - failures}/{len(expected)} 통과")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
