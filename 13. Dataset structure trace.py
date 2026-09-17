from ccdc.io import EntryReader  # 文档推荐的EntryReader初始化方式
import pandas as pd
import os
from tqdm import tqdm

INPUT_FILE = r"固体溶液refcode.xlsx"
OUTPUT_FILE = r"SS_info_result.xlsx"
STDERR_FILE = os.path.splitext(os.path.basename(__file__))[0] + '.stderr.txt'


def main():
    print("开始处理 CSD Refcode 信息检索（基于官方API）")
    print(f"输入文件: {INPUT_FILE}")
    print(f"输出文件: {OUTPUT_FILE}")

    # 按文档推荐方式初始化 EntryReader（替代之前的CrystalReader，更贴合publication属性调用）
    csd_reader = EntryReader('CSD')

    df = pd.read_excel(INPUT_FILE)
    refcodes = df.iloc[:, 0].tolist()
    years = []
    dois = []
    journal_names = []  # 期刊名
    article_titles = []  # 文章名（新增）
    missing_refcodes = []

    print("开始处理 Refcode 列表：")
    for idx, refcode in tqdm(enumerate(refcodes), total=len(refcodes), desc="处理进度"):
        row_num = idx + 2
        try:
            if pd.isna(refcode) or str(refcode).strip() == "":
                years.append("")
                dois.append("")
                journal_names.append("")
                article_titles.append("")
                continue

            # 按文档方式获取 Entry 对象
            entry = csd_reader.entry(str(refcode).strip().upper())
            cit = entry.publication  # 获取 Citation 类型的文献信息

            # 提取字段（严格遵循文档属性）
            year = str(cit.year) if cit.year else ""
            doi = cit.doi if cit.doi else ""
            journal_name = cit.journal.name if (hasattr(cit, 'journal') and cit.journal.name) else ""
            article_title = cit.title if (hasattr(cit, 'title') and cit.title) else ""  # 文章名

            # 存入结果
            years.append(year)
            dois.append(doi)
            journal_names.append(journal_name)
            article_titles.append(article_title)

        except Exception as e:
            error_msg = f"行号: {row_num}，Refcode: {refcode}，错误: {str(e)}"
            print(f"处理失败：{error_msg}")
            missing_refcodes.append(error_msg)
            years.append("")
            dois.append("")
            journal_names.append("")
            article_titles.append("")

    # 插入结果（可根据需求调整列顺序）
    df.insert(2, "Year", years)
    df.insert(3, "DOI", dois)
    df.insert(4, "Journal Name", journal_names)  # 期刊名
    df.insert(5, "Article Title", article_titles)  # 文章名

    # 保存文件
    df.to_excel(OUTPUT_FILE, index=False, engine='openpyxl')
    print(f"\n处理完成！结果已保存至：{OUTPUT_FILE}")

    # 错误日志
    if missing_refcodes:
        with open(STDERR_FILE, 'w', encoding='utf-8') as f:
            f.write("处理失败的Refcode列表：\n")
            for msg in missing_refcodes:
                f.write(f"{msg}\n")
        print(f"共 {len(missing_refcodes)} 条记录处理失败，详情见：{STDERR_FILE}")
    else:
        print("所有Refcode均处理成功，无错误记录")


if __name__ == '__main__':
    main()