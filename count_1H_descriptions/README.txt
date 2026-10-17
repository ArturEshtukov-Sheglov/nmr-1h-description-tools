README - count_1H_descriptions.py
=================================

Survey of published 1H NMR descriptions (Organic Letters, first half of 2026)
Version: 15 October 2026

Author: Artur V. Eshtukov-Shcheglov
The code was written with the assistance of Claude Opus 5.5 (Anthropic) and was
reviewed and tested by the author.


1. Purpose
----------
The script reads Supporting Information (SI) PDF files and lists every text
1H NMR description it finds, e.g.

    1H NMR (400 MHz, CDCl3) δ 7.41 (d, J = 6.9 Hz, 2H), ...

For each description it records the compound label written before it (e.g. 3a)
and flags descriptions that are probably not counted as new compounds:

    KNOWN?     the text near the description says "known compound",
               "in agreement with the literature", "matched the reported data", etc.
    CITATION   a literature citation (superscript number or [n]) follows the
               compound label or name
    REPEAT     the same label already has a description (several solvents,
               rotamers, repeated experiment); the compound is counted once

As a cross-check it also counts 13C NMR descriptions and HRMS "found" values.
HRMS data are usually given only for new compounds, so a large difference
between the number of unflagged 1H descriptions and the number of HRMS values
points to an SI that should be checked by hand.


2. IMPORTANT - the result must be checked by hand
--------------------------------------------------
The script is a screening aid, not a counting rule. In the survey, the number
of compounds (n) was decided by reading each SI according to the pass-1 rules
(compounds marked as known by the authors are excluded; a citation attached only
to a synthetic procedure does not make a compound known; peptides made by SPPS
are not counted; several descriptions of one compound count once).

Test against the manually determined counts for 13 articles:
    exact or within +/-2 compounds ..... 9 of 13 articles
    larger deviation ................... 4 of 13 articles (-4, +3, -11, +8)
The larger deviations came mainly from known compounds that were missed or
wrongly flagged. Treat the output as a list of candidates, never as the final
count.


3. Requirements
---------------
Python 3.9 or newer, and two packages:
    pymupdf   (reads the PDF files)
    openpyxl  (writes the Excel table)
Install them once:  python -m pip install -r requirements.txt  (from the main
folder; tested versions)  or  python -m pip install pymupdf openpyxl
If a package is missing, the script stops and says which one.


4. How to run
-------------
A) Windows (terminal in Visual Studio, VS Code or PowerShell)
   1. Put count_1H_descriptions.py in a folder.
   2. In the terminal, go to that folder:
          cd "C:\path\to\the\folder"
   3. Run the script with one or several SI files:
          python count_1H_descriptions.py "C:\path\to\SI_1.pdf" "C:\path\to\SI_2.pdf"
      (If the PDF is in the same folder, the file name alone is enough.
       A PDF can also be dragged from Explorer into the terminal window.)

B) Google Colab (nothing to install on your computer)
   1. Paste the whole script into a code cell and run it.
   2. Upload the SI PDF file(s) when asked.
   3. Download the result tables from the Files panel on the left.


5. Output
---------
On screen, for each SI file:
    number of 1H descriptions, 13C descriptions and HRMS "found" values;
    number of repeated labels and of descriptions flagged as possibly known;
    candidate n = descriptions without any flag (to be checked by hand);
    all labels in SI order ("?" = no label found before the description).

Two tables in the folder where the script was run (one row per 1H description):
    1H_descriptions.xlsx   for Excel (opens correctly in any language version)
    1H_descriptions.csv    standard CSV (comma-separated, UTF-8) for Python, R, etc.

Columns:
    file          SI file name
    index         number of the description in SI order
    label         compound label found before the description ("" = none)
    flags         KNOWN? / CITATION / REPEAT, empty if none
    description   full text of the description, from "1H NMR" up to the next
                  data block (13C/19F/31P/11B NMR, HRMS, IR, m.p., Rf, [α]);
                  at most 4000 characters

Both tables are overwritten every time the script is run. To keep the results
of several runs, rename the files, or give all PDFs in one command.


6. Troubleshooting
------------------
"Permission denied: 1H_descriptions..."
    The table is open in Excel. Close it and run the script again.

"python" or "pip" is not recognized
    Python is not installed or not on the PATH. Install Python from python.org
    and tick "Add Python to PATH" during installation.

Warning about "mkl-service" (from numpy)
    Harmless and unrelated to this script; the script hides it.
    Optional fix: pip install mkl-service

A description is cut or contains stray text
    PDF text extraction is not perfect (two-column layouts, figures inside the
    text, unusual fonts). Always compare with the PDF.

"?" instead of a label
    The label is written in a form the script does not recognise (e.g. only in
    bold, without parentheses). The description itself is still listed.


7. Limitations
--------------
- Descriptions given only in tables or as images are not found.
- Superscript citations are recognised only if the PDF stores them as
  superscript text; citations in other forms may be missed.
- The script does not judge whether a compound is new; it only marks candidates.
