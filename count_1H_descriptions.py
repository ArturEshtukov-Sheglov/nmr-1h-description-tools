# Counting 1H NMR descriptions in Supporting Information PDFs
# Survey of published 1H NMR descriptions (Organic Letters, H1 2026)
#
# Author:  Artur V. Eshtukov-Shcheglov
# Code written with the assistance of Claude Opus 5.5 (Anthropic);
# reviewed and tested by the author.
# Version: 15 October 2026
#
# IMPORTANT - the result must always be checked by hand.
#   Tested against the manually determined pass-1 counts for 13 articles:
#   the script agreed exactly or within +/-2 compounds in 9 of 13 articles,
#   and deviated by more in 4 articles (by -4, +3, -11 and +8), mainly because
#   "known" compounds were missed or wrongly flagged. Treat its output as a
#   list of candidates, never as the final count.
#
# What it does
#   For every SI PDF given, the script finds text 1H NMR descriptions
#   ("1H NMR (400 MHz, CDCl3) δ ..."), takes the compound label written just
#   before each description (e.g. "(3a)"), and flags descriptions that are
#   probably NOT counted under the pass-1 rules:
#     KNOWN   - "known compound", "matched the literature", "in agreement with
#               reported data" etc. near the description, or a literature
#               citation attached to the compound name/label;
#     REPEAT  - the same label already has a description (several solvents,
#               rotamers, a repeated scale-up) - the compound counts once.
#   It also counts 13C NMR descriptions and HRMS "found" values as a
#   cross-check: large differences between the three numbers point to papers
#   that should be checked by hand.
#
# What it does NOT do
#   The script is a screening aid. The final number of compounds n was
#   decided by reading the SI (pass-1 rules: compounds marked as known by the
#   authors are excluded; a citation attached only to a procedure does not
#   make a compound known; SPPS peptides are not counted; etc.). Every flag
#   must therefore be checked in the PDF.
#
# How to run in Google Colab
#   1) Run this cell. 2) Upload one or several SI PDFs when asked.
#   3) The script prints a summary per file and saves the table
#      (one row per 1H description, with the full text of the description)
#      as 1H_descriptions.xlsx and 1H_descriptions.csv.

import re, sys, csv, subprocess

try:
    import pymupdf as fitz  # PyMuPDF (pip install pymupdf)
except ImportError:
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", "pymupdf"], check=True)
    import pymupdf as fitz

# ---------------------------------------------------------------- text
def pdf_text(path):
    """Plain text of the whole PDF; page headers/footers removed; whitespace collapsed.
    Superscript text (literature citations such as "5:16") is kept as ^{16}."""
    doc = fitz.open(path)
    out = []
    for page in doc:
        H = page.rect.height
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines", []):
                y0, y1 = line["bbox"][1], line["bbox"][3]
                raw = "".join(sp["text"] for sp in line["spans"])
                if (y1 < 0.08 * H or y0 > 0.90 * H) and re.fullmatch(r"\s*-?\s?S?\s?\d{1,3}\s?-?\s*", raw):
                    continue                  # page number in the header or footer
                parts = []
                for span in line["spans"]:
                    t = span["text"]
                    if not t.strip():
                        parts.append(t); continue
                    if span["flags"] & 1 and re.fullmatch(r"\s*[\d,–\-\s]+\s*", t):
                        t = "^{" + t.strip() + "}"     # superscript digits = citation
                    parts.append(t)
                s = "".join(parts).strip()
                if not s or "Downloaded from" in s or re.fullmatch(r"-?S\s?\d+-?", s):
                    continue
                out.append(s)
    text = re.sub(r"\s+", " ", " ".join(out))
    # isotope superscripts ("^{1}H", "^{13}C", "^{19}F") are not citations
    text = re.sub(r"\^\{(\d{1,3})\}\s?(?=(?:H|C|F|P|B|N|Si|Sn|Pt|Se)\b|(?:H|C|F|P|B|N)[\s{-])", r"\1", text)
    return text

# ---------------------------------------------------------------- patterns
# 1H header: "1H NMR (", "1H-NMR (", "¹H NMR:", "1H NMR(400 MHz" ...
H1 = re.compile(r"(?<![0-9A-Za-z])(?:1\s?H|¹H)\s?-?\s?NMR\s?(?:\(|:|\[)")
C13 = re.compile(r"(?<![0-9A-Za-z])(?:13\s?C|¹³C)\s?(?:\{\s?(?:1\s?H|¹H)\s?\})?\s?-?\s?NMR\s?(?:\(|:|\[)")
HRMS_FOUND = re.compile(r"\bfound\b\s*:?\s*\d{2,4}\.\d{3,5}", re.I)

# a description follows the header: a delta sign and at least one "x.xx (" within ~200 characters
DESC_BODY = re.compile(r"δ.{0,200}?\d\.\d{1,3}\s?(?:–|-|—)?\s?(?:\d\.\d{1,3}\s?)?\(", re.S)

# compound label in parentheses, e.g. (3a), (5'), (S-2), (±)-12b is not taken
LABEL = re.compile(r"\(\s*((?:[A-Z]{1,2}-?)?\d{1,3}[a-z]{0,3}[′'']?(?:-[a-z]{1,2})?)\s*\)"
                   r"|\b[Cc]ompound\s+((?:[A-Z]{1,2}-?)?\d{1,3}[a-z]{0,3}[′'']?)")

KNOWN = re.compile(
    r"known compound|previously (?:reported|described|prepared|synthesi[sz]ed)"
    r"|(?:in (?:good )?(?:accordance|agreement)|consistent|identical) with (?:the |those |previously |literature|reported)"
    r"|match(?:ed|es|ing)? (?:well )?(?:with )?(?:those|the|literature|reported|previously)"
    r"|according to (?:the )?literature|literature data|reported data",
    re.I)
# citation right after a label or name: "(3a)12", "3a[12]", superscript digits are lost
# in text extraction, so "known" by citation can only be suggested, not decided
CITED_LABEL = re.compile(r"\)?\s?:?\s?(?:\^\{[\d,–\- ]+\}|\[\d{1,3}(?:[,–-]\d{1,3})*\])")
# a 1H description ends where the next data block begins
DESC_END = re.compile(r"(?<![0-9A-Za-z])(?:13\s?C|19\s?F|31\s?P|11\s?B|29\s?Si|1\s?H)\s?(?:\{[^}]*\})?\s?-?\s?NMR"
                      r"|\bHRMS\b|\bHR-?MS\b|\bIR\b\s?\(|\bFTIR\b|\bMS\s?\((?:ESI|EI|APCI)|\bm\.?p\.|\bR\s?f\b|\[α\]")

CAPTION = re.compile(r"(?:Figure|Fig\.|Spectrum|spectra)\s?S?\d*\.?\s.{0,40}$", re.I)

# ---------------------------------------------------------------- analysis
def full_description(S, start, limit=4000):
    """Text of the 1H description: from its header up to the next NMR header,
    HRMS, IR, m.p., Rf or optical rotation (at most `limit` characters)."""
    e = DESC_END.search(S, start + 8)
    end = e.start() if e and e.start() - start < limit else start + limit
    return S[start:end].strip()

def analyse(path):
    S = pdf_text(path)
    rows, seen = [], {}
    for m in H1.finditer(S):
        body = S[m.start(): m.start() + 400]
        if not DESC_BODY.search(body):
            continue                     # header without shifts (e.g. spectrum caption)
        before = S[max(0, m.start() - 60): m.start()]
        if CAPTION.search(before):
            continue                     # "Figure S5. 1H NMR (400 MHz) of 3a" - not a description
        # label: last "(xx)" in the 700 characters before the header
        window = S[max(0, m.start() - 700): m.start()]
        labs = list(LABEL.finditer(window))
        last = labs[-1] if labs else None
        label = (last.group(1) or last.group(2)) if last else ""
        ctx = S[max(0, m.start() - 400): m.start() + 1500]
        flags = []
        k = KNOWN.search(ctx)
        if k:
            flags.append("KNOWN?: " + k.group(0))
        if label:
            if CITED_LABEL.match(window, last.end()):
                flags.append("CITATION after label?")
            if label in seen:
                flags.append("REPEAT of description #%d" % seen[label])
            else:
                seen[label] = len(rows) + 1
        rows.append(dict(file=path.split("/")[-1], index=len(rows) + 1, label=label,
                         flags="; ".join(flags), description=full_description(S, m.start())))
    n13 = sum(1 for m in C13.finditer(S) if "δ" in S[m.start(): m.start() + 200])
    nhr = len(HRMS_FOUND.findall(S))
    return rows, n13, nhr

def summary(path, rows, n13, nhr):
    n1h = len(rows)
    labels = [r["label"] for r in rows]
    repeats = sum("REPEAT" in r["flags"] for r in rows)
    known = sum("KNOWN" in r["flags"] or "CITATION" in r["flags"] for r in rows)
    unflagged = sum(1 for r in rows if not r["flags"])
    print("=" * 70)
    print(path.split("/")[-1])
    print(f"  1H descriptions: {n1h} | 13C descriptions: {n13} | HRMS 'found': {nhr}")
    print(f"  repeated labels: {repeats} | flagged as possibly known: {known}")
    print(f"  candidate n (descriptions without any flag): {unflagged}  <- check by hand")
    print("  labels in SI order:", " ".join(l or "?" for l in labels))

# ---------------------------------------------------------------- run
if __name__ == "__main__":
    files = sys.argv[1:]
    if not files:
        try:
            from google.colab import files as colab_files   # Colab: upload dialog
            files = list(colab_files.upload().keys())
        except ImportError:
            sys.exit("Usage: python count_1H_descriptions.py SI1.pdf [SI2.pdf ...]")
    all_rows = []
    for f in files:
        rows, n13, nhr = analyse(f)
        summary(f, rows, n13, nhr)
        all_rows += rows
    fields = ["file", "index", "label", "flags", "description"]
    # 1) Excel table - opens correctly in any language version of Excel
    import warnings   # openpyxl imports numpy if installed; hide numpy's unrelated MKL warning
    warnings.filterwarnings("ignore", message=".*mkl-service.*")
    try:
        import openpyxl
    except ImportError:
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", "openpyxl"], check=True)
        import openpyxl
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = "1H descriptions"
    ws.append(fields)
    for r in all_rows:
        ws.append([r[f] for f in fields])
    for col, width in zip("ABCDE", (24, 7, 10, 40, 150)):
        ws.column_dimensions[col].width = width
    for row in ws.iter_rows(min_row=2, min_col=5, max_col=5):
        row[0].alignment = openpyxl.styles.Alignment(wrap_text=True, vertical="top")
    try:
        wb.save("1H_descriptions.xlsx")
        # 2) standard CSV (comma-separated, UTF-8) for Python, R and other programs
        with open("1H_descriptions.csv", "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader(); w.writerows(all_rows)
    except PermissionError:
        sys.exit("\nCannot write the table: 1H_descriptions.xlsx or .csv is open in another "
                 "program (e.g. Excel). Close it and run the script again.")
    print("\nTables saved: 1H_descriptions.xlsx (for Excel) and 1H_descriptions.csv"
          " (one row per 1H description)")
