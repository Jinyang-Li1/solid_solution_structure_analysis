import pandas as pd
from crossref.restful import Works
from tqdm import tqdm
import time


def get_article_title(doi, max_retries=3):
    if not doi or pd.isna(doi):
        return None

    works = Works()
    for retry in range(max_retries):
        try:
            result = works.doi(doi)
            if result and 'title' in result and len(result['title']) > 0:
                return result['title'][0]
            else:
                return "Title not found"
        except Exception as e:
            if retry < max_retries - 1:
                time.sleep(1)
                continue
            return f"Query failed: {str(e)[:50]}..."


def main():
    input_file = "SS_info_result.xlsx"
    output_file = "SS_info_result_with_titles.xlsx"

    try:
        print(f"Reading file: {input_file}")
        df = pd.read_excel(input_file, sheet_name='Sheet1')

        required_columns = ['DOI', 'Article Title']
        for col in required_columns:
            if col not in df.columns:
                raise ValueError(f"Missing required column: {col}")

        print(f"Found {len(df)} records, starting query...\n")

        for idx, row in tqdm(enumerate(df.itertuples(index=False)), total=len(df), desc="Progress"):
            doi = row.DOI
            title = get_article_title(doi)
            df.at[idx, 'Article Title'] = title
            print(
                f"Row {idx + 1}/{len(df)} | DOI: {doi if pd.notna(doi) else 'No DOI'} | Title: {title if title else 'None'}")

        df.to_excel(output_file, index=False)
        print(f"\nCompleted! Results saved to: {output_file}")

    except Exception as e:
        print(f"Error occurred: {str(e)}")


if __name__ == "__main__":
    print("Install dependencies before first run:")
    print("pip install pandas crossrefapi tqdm openpyxl")
    print("----------------------------------------")
    main()