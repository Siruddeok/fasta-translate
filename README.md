# fasta-translate

FASTA 파일의 DNA/RNA 서열을 검증하고 아미노산으로 번역해 CSV로 저장하는 Biopython 스크립트입니다.
실제 공개 indel 돌연변이(ClinVar/문헌, NCBI RefSeq 기반)로 동작을 검증했습니다.

## 설치

```bash
pip install -r requirements.txt
```

## 사용법

```bash
python fasta_translate.py input.fasta -o result.csv
```

| 옵션 | 설명 |
|---|---|
| `-o`, `--output` | 출력 CSV 경로 (기본: `<입력파일명>_translated.csv`) |
| `--start frame1` | 1번 염기부터 서열 끝까지 번역 (기본값, CDS 입력용) |
| `--start first-atg` | 첫 ATG부터 그 틀의 첫 종결 코돈까지 번역 (5'UTR이 있는 mRNA 입력용) |
| `--start longest-orf` | 정방향 3개 틀에서 가장 긴 ORF를 번역 (정체를 모르는 서열용) |
| `--alt-starts` | `first-atg`/`longest-orf`에서 유전 코드의 대체 개시 코돈(GTG, TTG 등)도 찾음 |
| `--table N` | NCBI 유전 코드 번호 (기본 1: 표준, 2: 척추동물 미토콘드리아, 11: 세균) |
| `--strict` | 경고도 오류로 취급해 완전한 CDS만 유효로 인정 |

## 검증 규칙

번역이 불가능한 경우는 **오류**(무효, 번역하지 않음), 번역은 가능하지만 확인이 필요한 경우는 **경고**로 기록합니다.

- 오류: 빈 서열, 허용되지 않은 문자, T와 U 혼재, 3염기 미만, 첫 코돈이 종결 코돈, 개시 코돈을 찾지 못함
- 경고: 길이가 3의 배수가 아님, 개시 코돈 없음 / 대체 개시 코돈으로 시작, 종결 코돈으로 끝나지 않음,
  중간 종결 코돈, 모호 염기 N, 갭 문자 제거, 중복 ID

RNA는 U→T로 바꿔 처리하며, 개시·종결 코돈 목록은 `--table`의 유전 코드를 따릅니다.

## 출력 CSV 열

`ID, Description, Seq_Type, Length, GC_Content, Start_Pos, End_Pos, Is_Valid, Errors, Warnings,
Protein_Full, Protein_To_Stop, Protein_Length`

- `Protein_Full`: 번역 구간 전체 번역 (`*` = 종결 코돈)
- `Protein_To_Stop`: 첫 종결 코돈 직전까지의 단백질
- `Start_Pos`/`End_Pos`: 번역 구간 (1-based, 갭 제거 후 좌표)

## 검증 (`validation/`)

### 1. 실제 indel 돌연변이

NCBI RefSeq에서 CDS를 내려받아 HGVS 표기대로 indel을 적용하고, HGVS p. 표기에서 계산한 단백질 길이와 비교합니다.

| 유전자 | 변이 | 질환 | 종류 |
|---|---|---|---|
| HBB | c.126_129del | 베타지중해빈혈 | −4 nt, frameshift |
| CFTR | c.1521_1523del (ΔF508) | 낭포성 섬유증 | −3 nt, in-frame |
| BRCA1 | c.68_69del (185delAG) | 유전성 유방암 | −2 nt, frameshift |
| BRCA2 | c.5946del (6174delT) | 유전성 유방암 | −1 nt, frameshift |
| HEXA | c.1274_1277dup | 테이-삭스병 | +4 nt, frameshift |
| CCR5 | c.554_585del (Δ32) | HIV 저항성 | −32 nt, frameshift |

```bash
python validation/build_indel_fasta.py      # RefSeq 다운로드 (data/에 캐시) + FASTA 생성
python fasta_translate.py validation/indel_test.fasta -o validation/indel_translated.csv
python validation/compare_results.py        # 결과: 12/12 통과
```

### 2. 개시 코돈으로 시작하지 않는 서열

5'UTR 포함 mRNA, 틀이 어긋난 CDS, CDS 내부 조각, 첫 코돈이 종결 코돈, GTG 시작 등 8가지 경우.

```bash
python validation/build_no_start_fasta.py
python fasta_translate.py validation/no_start_test.fasta -o validation/ns_first-atg_t1.csv --start first-atg
```

`validation/ns_<모드>_t<유전코드>.csv`에 모드별 결과가 있습니다.

## 한계

- 기준 서열과 비교하지 않으므로, ΔF508처럼 3의 배수 길이의 in-frame indel은 경고 없이 정상 서열로 처리됩니다.
- `longest-orf`는 정방향 가닥만 검색합니다.
