import pandas as pd
from crossref.restful import Works
from tqdm import tqdm  # 用于显示进度条
import time


def get_article_title(doi, max_retries=3):
    """根据DOI获取文章标题，包含重试机制"""
    if not doi or pd.isna(doi):
        return None

    works = Works()
    for retry in range(max_retries):
        try:
            result = works.doi(doi)
            if result and 'title' in result and len(result['title']) > 0:
                return result['title'][0]
            else:
                return f"未找到标题"
        except Exception as e:
            if retry < max_retries - 1:
                time.sleep(1)
                continue
            return f"查询失败: {str(e)[:50]}..."  # 简化错误信息长度


def main():
    # 本地文件路径（请替换为你的实际路径）
    input_file = "SS_info_result.xlsx"
    output_file = "SS_info_result_with_titles.xlsx"

    try:
        print(f"正在读取文件: {input_file}")
        df = pd.read_excel(input_file, sheet_name='Sheet1')

        required_columns = ['DOI', 'Article Title']
        for col in required_columns:
            if col not in df.columns:
                raise ValueError(f"缺少必要的列: {col}")

        print(f"共发现 {len(df)} 条记录，开始查询并打印每行结果...\n")

        # 逐行处理并打印结果
        for idx, row in tqdm(enumerate(df.itertuples(index=False)), total=len(df), desc="处理进度"):
            doi = row.DOI
            # 获取标题
            title = get_article_title(doi)
            # 保存标题到DataFrame
            df.at[idx, 'Article Title'] = title
            # 打印当前行结果（索引从1开始更直观）
            print(
                f"行 {idx + 1}/{len(df)} | DOI: {doi if pd.notna(doi) else '无DOI'} | 标题: {title if title else '无'}")

        # 保存结果
        df.to_excel(output_file, index=False)
        print(f"\n处理完成！结果已保存到: {output_file}")

    except Exception as e:
        print(f"发生错误: {str(e)}")


if __name__ == "__main__":
    print("首次运行请先安装依赖库：")
    print("pip install pandas crossrefapi tqdm openpyxl")
    print("----------------------------------------")
    main()