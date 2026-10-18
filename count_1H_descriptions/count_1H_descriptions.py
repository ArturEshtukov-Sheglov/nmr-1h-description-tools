"""Count text 1H NMR descriptions in Supporting Information (SI) PDF files.

Lists every description such as "1H NMR (400 MHz, CDCl3) δ 7.41 (d, J = 6.9 Hz, 2H), ...",
the compound label written before it, and flags probably known (KNOWN?, CITATION) or
repeated (REPEAT) compounds. A screening aid: every flag must be checked in the PDF.
Purpose, test results, output columns and limitations: README.txt in this folder.

    python count_1H_descriptions.py SI_1.pdf [SI_2.pdf ...] [-o 1H_descriptions]

This file is self-contained on purpose: it can be pasted into one Google Colab cell.
extract_names.py and check_descriptions.py import the description patterns from it.

Author: Artur V. Eshtukov-Shcheglov. Code written with the assistance of Claude Opus 5.5
(Anthropic); reviewed and tested by the author. Version 1.1 (October 2026).
"""

import csv
import os
import re
import sys

try:
    import pymupdf as fitz  # PyMuPDF
except ImportError:
    sys.exit("The package pymupdf is missing. Install the packages with:\n    python -m pip install -r requirements.txt\n(or: python -m pip install pymupdf openpyxl rdkit)")


# ===== Constants =====

# Page layout: a line in the top 8 % or the bottom 10 % of a page that is only a page number
# ("7", "S12", "- 3 -") is a page header or footer and is dropped.
PAGE_HEADER_ZONE = 0.08          # fraction of the page height from the top
PAGE_FOOTER_ZONE = 0.90          # lines starting below this fraction of the height
PAGE_NUMBER_RE = (
    r"\s*-?\s?"                  # optional dash before the number: "- 3 -"
    r"S?\s?"                     # optional "S" of SI page numbers: "S12"
    r"\d{1,3}"                   # the page number itself (1-999)
    r"\s?-?\s*"                  # optional dash after the number
)
# a whole line that is an SI page number ("S12", "-S3-"), wherever it is on the page
SI_PAGE_LINE_RE = r"-?S\s?\d+-?"
# superscript digits (PyMuPDF span flag 1) such as "5", "12,13", "3-5": a literature citation
SUPERSCRIPT_FLAG = 1
CITATION_DIGITS_RE = r"\s*[\d,–\-\s]+\s*"
# isotope numbers written as superscripts ("^{1}H", "^{13}C") are not citations: written
# back as plain "1H", "13C"; a number of 1-3 digits followed by a nucleus symbol
ISOTOPE_SUPERSCRIPT_RE = (
    r"\^\{(\d{1,3})\}\s?"
    r"(?=(?:H|C|F|P|B|N|Si|Sn|Pt|Se)\b"   # followed by a whole element symbol
    r"|(?:H|C|F|P|B|N)[\s{-])"            # or by H/C/F/P/B/N and a space, "{" or "-" ("1H-NMR")
)
# a δ typed in the Symbol font is extracted as the private-use character U+F064
SYMBOL_FONT_DELTA = ""

# Text windows around a 1H header, in characters of the extracted text
BODY_WINDOW = 400                # after the header: must contain the shifts of a description
CAPTION_WINDOW = 60              # before the header: "Figure S5." here means a spectrum caption
LABEL_WINDOW = 700               # before the header: the compound label is looked for here
KNOWN_WINDOW_BEFORE = 400        # before the header: "known compound" etc. is looked for here
KNOWN_WINDOW_AFTER = 1500        # after the header: the same
MAX_DESCRIPTION = 4000           # longest description kept
DESC_END_SKIP = 8                # the header itself ("1H NMR (") is not a data-block end
C13_BODY_WINDOW = 200            # after a 13C header: must contain a δ to count as a description


# ===== Patterns =====

# 1H header: "1H NMR (", "1H-NMR (", "¹H NMR:", "1H NMR(400 MHz", "1H NMR ["
H1 = re.compile(
    r"(?<![0-9A-Za-z])"          # not inside a longer word or number ("21H", "D1H")
    r"(?:1\s?H|¹H)"              # "1H", "1 H" or "¹H"
    r"\s?-?\s?NMR\s?"            # "NMR", "-NMR", " NMR"
    r"(?:\(|:|\[)"               # followed by "(", ":" or "["
)
# 13C header, also "13C{1H} NMR (": counted as a cross-check only
C13 = re.compile(
    r"(?<![0-9A-Za-z])"
    r"(?:13\s?C|¹³C)"
    r"\s?(?:\{\s?(?:1\s?H|¹H)\s?\})?"     # optional "{1H}" (proton-decoupled)
    r"\s?-?\s?NMR\s?"
    r"(?:\(|:|\[)"
)
# HRMS "found 312.1234": counted as a cross-check only
HRMS_FOUND = re.compile(r"\bfound\b\s*:?\s*\d{2,4}\.\d{3,5}", re.I)

# The header is followed by a real description (not by a spectrum caption) if one of two forms
# is found within BODY_WINDOW characters:
DESC_BODY = re.compile(
    # (1) a δ and, within 200 characters, a shift or shift range followed by "(":
    #     "δ 7.41 (d, ...", "δ 7.35-7.28 (m, ..."
    r"δ.{0,200}?\d\.\d{1,3}\s?(?:–|-|—)?\s?(?:\d\.\d{1,3}\s?)?\("
    # (2) without δ: "1H NMR (300 MHz, CDCl3): 7.35-7.28 (m, 2H)"
    r"|^[^)\]]{0,80}[)\]]?"                        # the rest of the header, up to its ")" or "]"
    r"\s?:?\s?"                                    # optional ":"
    r"(?:[A-Za-z][^()]{0,60}?[:.]\s?)?"            # optional word(s) such as "major isomer:"
    r"\d{1,2}\.\d{1,3}\s?"                         # the first shift
    r"(?:[–\-—−]\s?\d{1,2}\.\d{1,3}\s?)?"          # optional end of a range
    r"\(\s?(?:app(?:arent)?\.?\s|br\.?\s|ps\s)?"   # "(" and optional "app", "br", "ps"
    r"(?:s|d|t|q|m|p|h|br|sept|hept|quint|pent|sext|[dtq]{2,5})\b",   # a multiplicity
    re.S)

# Compound label in parentheses before the description: (3a), (5'), (S-2), (12b-ii);
# or "compound 3a". "(±)-12b" is not taken.
LABEL = re.compile(
    r"\(\s*("
    r"(?:[A-Z]{1,2}-?)?"         # optional letter prefix: "S-", "SI", "A"
    r"\d{1,3}"                   # number
    r"[a-z]{0,3}"                # optional letters: "a", "ab"
    r"[′'’]?"                    # optional prime
    r"(?:-[a-z]{1,2})?"          # optional suffix: "-ii"
    r")\s*\)"
    r"|\b[Cc]ompound\s+((?:[A-Z]{1,2}-?)?\d{1,3}[a-z]{0,3}[′'’]?)"
)

# Text near the description saying that the compound is known
KNOWN = re.compile(
    r"known compound|previously (?:reported|described|prepared|synthesi[sz]ed)"
    r"|(?:in (?:good )?(?:accordance|agreement)|consistent|identical) with (?:the |those |previously |literature|reported)"
    r"|match(?:ed|es|ing)? (?:well )?(?:with )?(?:those|the|literature|reported|previously)"
    r"|according to (?:the )?literature|literature data|reported data",
    re.I)
# A citation right after a label or name: "(3a)^{12}", "3a[12]", "3a [12-14]". Superscript digits
# can be lost in text extraction, so a citation can only suggest that a compound is known.
CITED_LABEL = re.compile(
    r"\)?\s?:?\s?"
    r"(?:\^\{[\d,–\- ]+\}"                         # superscript citation: ^{12} or ^{12,13}
    r"|\[\d{1,3}(?:[,–-]\d{1,3})*\])"              # citation in brackets: [12], [12,13], [12-14]
)
# A 1H description ends where the next data block begins: another NMR header, HRMS, IR, MS,
# melting point, Rf or optical rotation
DESC_END = re.compile(
    r"(?<![0-9A-Za-z])(?:13\s?C|19\s?F|31\s?P|11\s?B|29\s?Si|1\s?H)\s?(?:\{[^}]*\})?\s?-?\s?NMR"
    r"|\bHRMS\b|\bHR-?MS\b"
    r"|\bIR\b\s?\(|\bFTIR\b"
    r"|\bMS\s?\((?:ESI|EI|APCI)"
    r"|\bm\.?p\.|\bR\s?f\b|\[α\]"
)
# Text just before a header that makes it a spectrum caption: "Figure S5. 1H NMR (400 MHz) of 3a"
CAPTION = re.compile(r"(?:Figure|Fig\.|Spectrum|spectra)\s?S?\d*\.?\s.{0,40}$", re.I)

OUTPUT_FIELDS = ["file", "index", "label", "flags", "description"]
COLUMN_WIDTHS = (24, 7, 10, 40, 150)
MISSING_OPENPYXL = ("The package openpyxl is missing. Install the packages with:\n    python -m pip install -r requirements.txt\n"
                    "(or: python -m pip install pymupdf openpyxl rdkit)")


# ===== Reading the PDF =====

def is_page_number(line, raw, page_height):
    """True if a text line is a page number in the page header or footer."""
    y0, y1 = line["bbox"][1], line["bbox"][3]
    in_margin = y1 < PAGE_HEADER_ZONE * page_height or y0 > PAGE_FOOTER_ZONE * page_height
    return in_margin and re.fullmatch(PAGE_NUMBER_RE, raw)


def line_text(line):
    """Text of one line; superscript citation digits are written as ^{16}."""
    parts = []
    for span in line["spans"]:
        t = span["text"]
        if not t.strip():
            parts.append(t)
            continue
        if span["flags"] & SUPERSCRIPT_FLAG and re.fullmatch(CITATION_DIGITS_RE, t):
            t = "^{" + t.strip() + "}"
        parts.append(t)
    return "".join(parts).strip()


def pdf_text(path):
    """Plain text of the whole PDF; page headers/footers removed; whitespace collapsed.
    Superscript text (literature citations such as "5:16") is kept as ^{16}."""
    doc = fitz.open(path)
    out = []
    for page in doc:
        page_height = page.rect.height
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines", []):
                raw = "".join(sp["text"] for sp in line["spans"])
                if is_page_number(line, raw, page_height):
                    continue

                s = line_text(line)
                if not s or "Downloaded from" in s or re.fullmatch(SI_PAGE_LINE_RE, s):
                    continue
                out.append(s)

    text = re.sub(r"\s+", " ", " ".join(out))
    text = text.replace(SYMBOL_FONT_DELTA, "δ")
    text = re.sub(ISOTOPE_SUPERSCRIPT_RE, r"\1", text)
    return text


# ===== Finding the descriptions =====

def is_description_header(S, pos):
    """A 1H header at pos starts a description: shifts follow, and it is not a spectrum caption."""
    return (bool(DESC_BODY.search(S[pos: pos + BODY_WINDOW]))
            and not CAPTION.search(S[max(0, pos - CAPTION_WINDOW): pos]))


def full_description(S, start, limit=MAX_DESCRIPTION):
    """Text of the 1H description: from its header up to the next NMR header,
    HRMS, IR, m.p., Rf or optical rotation (at most `limit` characters)."""
    e = DESC_END.search(S, start + DESC_END_SKIP)
    end = e.start() if e and e.start() - start < limit else start + limit
    return S[start:end].strip()


def label_before(S, pos):
    """The last label "(xx)" or "compound xx" in the LABEL_WINDOW characters before pos.
    Returns (label, match, window); label is "" if none is found."""
    window = S[max(0, pos - LABEL_WINDOW): pos]
    labs = list(LABEL.finditer(window))
    last = labs[-1] if labs else None
    label = (last.group(1) or last.group(2)) if last else ""
    return label, last, window


def analyse(path):
    """All 1H descriptions of one SI: (rows, number of 13C descriptions, number of HRMS values)."""
    S = pdf_text(path)
    rows, seen = [], {}
    for m in H1.finditer(S):
        if not is_description_header(S, m.start()):
            continue                     # no shifts follow, or a spectrum caption

        label, last, window = label_before(S, m.start())
        context = S[max(0, m.start() - KNOWN_WINDOW_BEFORE): m.start() + KNOWN_WINDOW_AFTER]

        flags = []
        k = KNOWN.search(context)
        if k:
            flags.append("KNOWN?: " + k.group(0))
        if label:
            if CITED_LABEL.match(window, last.end()):
                flags.append("CITATION after label?")
            if label in seen:
                flags.append("REPEAT of description #%d" % seen[label])
            else:
                seen[label] = len(rows) + 1

        rows.append(dict(file=os.path.basename(path.replace("\\", "/")), index=len(rows) + 1, label=label,
                         flags="; ".join(flags), description=full_description(S, m.start())))

    n13 = sum(1 for m in C13.finditer(S) if "δ" in S[m.start(): m.start() + C13_BODY_WINDOW])
    nhr = len(HRMS_FOUND.findall(S))
    return rows, n13, nhr


# ===== Output =====

def summary(path, rows, n13, nhr):
    """Print the counts of one SI file."""
    n1h = len(rows)
    labels = [r["label"] for r in rows]
    repeats = sum("REPEAT" in r["flags"] for r in rows)
    known = sum("KNOWN" in r["flags"] or "CITATION" in r["flags"] for r in rows)
    unflagged = sum(1 for r in rows if not r["flags"])

    print("=" * 70)
    print(os.path.basename(path.replace("\\", "/")))
    print(f"  1H descriptions: {n1h} | 13C descriptions: {n13} | HRMS 'found': {nhr}")
    print(f"  repeated labels: {repeats} | flagged as possibly known: {known}")
    print(f"  candidate n (descriptions without any flag): {unflagged}  <- check by hand")
    print("  labels in SI order:", " ".join(l or "?" for l in labels))


def save_tables(all_rows, base):
    """Write base.xlsx (for Excel) and base.csv (UTF-8, for other programs)."""
    import warnings   # openpyxl imports numpy if installed; hide numpy's unrelated MKL warning
    warnings.filterwarnings("ignore", message=".*mkl-service.*")
    try:
        import openpyxl
    except ImportError:
        sys.exit(MISSING_OPENPYXL)

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "1H descriptions"
    ws.append(OUTPUT_FIELDS)
    for r in all_rows:
        ws.append([r[f] for f in OUTPUT_FIELDS])
    for col, width in zip("ABCDE", COLUMN_WIDTHS):
        ws.column_dimensions[col].width = width
    for row in ws.iter_rows(min_row=2, min_col=5, max_col=5):
        row[0].alignment = openpyxl.styles.Alignment(wrap_text=True, vertical="top")

    try:
        wb.save(base + ".xlsx")
        with open(base + ".csv", "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=OUTPUT_FIELDS)
            w.writeheader()
            w.writerows(all_rows)
    except PermissionError:
        sys.exit(f"\nCannot write the table: {base}.xlsx or .csv is open in another "
                 "program (e.g. Excel). Close it and run the script again.")
    print(f"\nTables saved: {base}.xlsx (for Excel) and {base}.csv"
          " (one row per 1H description)")


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Find all text 1H NMR descriptions in SI PDF files, with "
                                 "compound labels and flags for probably known or repeated compounds.")
    ap.add_argument("pdf", nargs="*", help="SI file(s) in PDF format (in Google Colab: none - an upload "
                                           "dialog opens)")
    ap.add_argument("-o", dest="out", default="1H_descriptions",
                    help="name of the output tables without extension (default: 1H_descriptions)")
    # in a Jupyter/Colab cell sys.argv holds the kernel's own arguments ("-f kernel.json"): ignore them
    a = ap.parse_args([] if "ipykernel" in sys.modules else None)
    files, base = a.pdf, re.sub(r"\.(xlsx|csv)$", "", a.out, flags=re.I)

    if not files:
        try:
            from google.colab import files as colab_files   # Colab: upload dialog
            files = list(colab_files.upload().keys())
        except ImportError:
            ap.error("give at least one PDF file")

    all_rows = []
    for f in files:
        rows, n13, nhr = analyse(f)
        summary(f, rows, n13, nhr)
        all_rows += rows
    save_tables(all_rows, base)


if __name__ == "__main__":
    main()
