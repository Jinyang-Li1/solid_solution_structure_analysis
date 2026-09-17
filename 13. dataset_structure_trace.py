from ccdc.io import EntryReader
import pandas as pd
import os
from tqdm import tqdm

INPUT_FILE = r"solid_solution_refcode.xlsx"
OUTPUT_FILE = r"SS_info_result.xlsx"
STDERR_FILE = os.path.splitext(os.path.basename(__file__))[0] + '.stderr.txt'


def main():
    print("Starting CSD Refcode metadata retrieval via CCDC API")
    print(f"Input file: {INPUT_FILE}")
    print(f"Output file: {OUTPUT_FILE}")

    csd_reader = EntryReader('CSD')

    df = pd.read_excel(INPUT_FILE)
    refcodes = df.iloc[:, 0].tolist()
    years = []
    dois = []
    journal_names = []
    article_titles = []
    missing_refcodes = []

    print("Processing Refcode list:")
    for idx, refcode in tqdm(enumerate(refcodes), total=len(refcodes), desc="Progress"):
        row_num = idx + 2
        try:
            if pd.isna(refcode) or str(refcode).strip() == "":
                years.append("")
                dois.append("")
                journal_names.append("")
                article_titles.append("")
                continue

            entry = csd_reader.entry(str(refcode).strip().upper())
            cit = entry.publication

            year = str(cit.year) if cit.year else ""
            doi = cit.doi if cit.doi else ""
            journal_name = cit.journal.name if (hasattr(cit, 'journal') and cit.journal.name) else ""
            article_title = cit.title if (hasattr(cit, 'title') and cit.title) else ""

            years.append(year)
            dois.append(doi)
            journal_names.append(journal_name)
            article_titles.append(article_title)

        except Exception as e:
            error_msg = f"Row: {row_num}, Refcode: {refcode}, Error: {str(e)}"
            print(f"Failed: {error_msg}")
            missing_refcodes.append(error_msg)
            years.append("")
            dois.append("")
            journal_names.append("")
            article_titles.append("")

    df.insert(2, "Year", years)
    df.insert(3, "DOI", dois)
    df.insert(4, "Journal Name", journal_names)
    df.insert(5, "Article Title", article_titles)

    df.to_excel(OUTPUT_FILE, index=False, engine='openpyxl')
    print(f"\nCompleted! Results saved to: {OUTPUT_FILE}")

    if missing_refcodes:
        with open(STDERR_FILE, 'w', encoding='utf-8') as f:
            f.write("Failed Refcode records:\n")
            for msg in missing_refcodes:
                f.write(f"{msg}\n")
        print(f"{len(missing_refcodes)} records failed. See log: {STDERR_FILE}")
    else:
        print("All Refcodes processed successfully, no errors")


if __name__ == '__main__':
    main()