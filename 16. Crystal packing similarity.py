import pandas as pd
import os
import csv
import sys
from collections import defaultdict

import ccdc.io
import ccdc.crystal


class RealTimeReporter:
    """实时写入包含晶胞参数的详细结果"""

    def __init__(self, detailed_file_path):
        self.detailed_file = detailed_file_path
        self.success_pairs = defaultdict(bool)  # 行级汇总用
        self._init_detailed_file()  # 初始化文件并写表头

    def _init_detailed_file(self):
        """初始化详细文件，写入包含晶胞参数的表头"""
        with open(self.detailed_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow([
                # 基础信息
                "group_id", "filename1", "filename2",
                "refcode1", "refcode2",
                # refcode1的晶胞参数
                "spacegroup1", "a1", "b1", "c1", "alpha1", "beta1", "gamma1",
                # refcode2的晶胞参数
                "spacegroup2", "a2", "b2", "c2", "alpha2", "beta2", "gamma2",
                # 比对结果
                "匹配分子数量", "RMSD", "状态", "备注（失败原因）"
            ])

    def _get_cell_params(self, entry):
        """提取单个entry的晶胞参数，返回元组"""
        try:
            crystal = entry.crystal
            spacegroup = crystal.spacegroup_number_and_setting[0]  # 空间群编号
            a = crystal.cell_lengths.a
            b = crystal.cell_lengths.b
            c = crystal.cell_lengths.c
            alpha = crystal.cell_angles.alpha
            beta = crystal.cell_angles.beta
            gamma = crystal.cell_angles.gamma
            return (spacegroup, round(a, 4), round(b, 4), round(c, 4),
                    round(alpha, 2), round(beta, 2), round(gamma, 2))
        except Exception as e:
            return ("获取失败", "", "", "", "", "", "")  # 提取失败时返回空

    def add_success(self, group_id, filename1, filename2, entry1, entry2, r1, r2, nmatch, rmsd):
        """实时写入成功比对的结果（含晶胞参数）"""
        # 提取晶胞参数
        cell1 = self._get_cell_params(entry1)  # (spacegroup1, a1, b1, c1, alpha1, beta1, gamma1)
        cell2 = self._get_cell_params(entry2)  # (spacegroup2, a2, b2, c2, alpha2, beta2, gamma2)

        # 判定状态
        status = "同构（15/15）" if nmatch == 15 else f"不同构（{nmatch}/15）"
        if nmatch == 15:
            self.success_pairs[group_id] = True

        # 实时追加写入
        with open(self.detailed_file, 'a', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow([
                group_id, filename1, filename2,
                r1, r2,
                # refcode1的晶胞参数
                cell1[0], cell1[1], cell1[2], cell1[3], cell1[4], cell1[5], cell1[6],
                # refcode2的晶胞参数
                cell2[0], cell2[1], cell2[2], cell2[3], cell2[4], cell2[5], cell2[6],
                # 比对结果
                nmatch, round(rmsd, 4), status, ""
            ])
            f.flush()

    def add_failure(self, group_id, filename1, filename2, r1, r2, reason, entry1=None, entry2=None):
        """实时写入失败比对的结果（尽可能提取晶胞参数）"""
        # 尝试提取可用的晶胞参数（若entry存在）
        cell1 = ("", "", "", "", "", "", "")
        cell2 = ("", "", "", "", "", "", "")
        if entry1:
            cell1 = self._get_cell_params(entry1)
        if entry2:
            cell2 = self._get_cell_params(entry2)

        with open(self.detailed_file, 'a', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow([
                group_id, filename1, filename2,
                r1, r2,
                # refcode1的晶胞参数（可能为空）
                cell1[0], cell1[1], cell1[2], cell1[3], cell1[4], cell1[5], cell1[6],
                # refcode2的晶胞参数（可能为空）
                cell2[0], cell2[1], cell2[2], cell2[3], cell2[4], cell2[5], cell2[6],
                # 比对结果
                "", "", "失败", reason
            ])
            f.flush()

    def get_summary(self, all_group_info):
        """生成行级汇总结果（包含identifier1和identifier2）"""
        summary = []
        # all_group_info格式：group_id -> (filename1, identifier1, filename2, identifier2)
        for group_id in sorted(all_group_info.keys()):
            filename1, identifier1, filename2, identifier2 = all_group_info[group_id]
            is_iso = "是" if self.success_pairs.get(group_id, False) else "否"
            summary.append((group_id, filename1, identifier1, filename2, identifier2, is_iso))
        return summary


class PairAnalysis:
    def __init__(self, pairs_with_meta, reporter):
        self.pairs_with_meta = pairs_with_meta  # (group_id, filename1, filename2, r1, r2)
        self.reporter = reporter
        if not self.pairs_with_meta:
            raise ValueError("未找到需要比对的refcode对")

    def similarity(self, ref_crystal, query_crystal):
        try:
            ps = ccdc.crystal.PackingSimilarity()
            ps.settings.ignore_hydrogen_positions = True
            ps.settings.ignore_hydrogen_counts = False
            ps.settings.ignore_bond_types = False
            ps.settings.ignore_bond_counts = True
            ps.settings.match_entire_packing_shell = False
            ps.settings.allow_molecular_differences = True
            ps.settings.ignore_smallest_components = True
            ps.settings.packing_shell_size = 15
            ps.settings.timeout_ms = 20000  # 延长超时时间到20秒，减少超时报错
            result = ps.compare(ref_crystal, query_crystal)
            return (result.nmatched_molecules, result.rmsd) if result else (None, "比对超时")
        except Exception as e:
            return (None, f"计算异常: {str(e)}")

    def run(self):
        reader = ccdc.io.EntryReader('CSD')
        total_pairs = len(self.pairs_with_meta)
        for idx, (group_id, filename1, filename2, r1, r2) in enumerate(self.pairs_with_meta):
            # 控制台实时进度
            progress = (idx + 1) / total_pairs * 100
            print(f"[进度: {progress:.1f}%] 处理第 {idx + 1}/{total_pairs} 对: group_id={group_id}, {r1} vs {r2}",
                  end='\r')

            entry1 = None
            entry2 = None
            try:
                # 读取晶体结构（保留entry对象用于提取晶胞参数）
                entry1 = reader.entry(r1)
                entry2 = reader.entry(r2)

                # 计算比对结果
                nmatch, rmsd_or_reason = self.similarity(entry1.crystal, entry2.crystal)
                if nmatch is None:
                    # 比对失败（但entry可能有效，尝试提取晶胞参数）
                    self.reporter.add_failure(
                        group_id, filename1, filename2, r1, r2,
                        reason=rmsd_or_reason, entry1=entry1, entry2=entry2
                    )
                else:
                    # 比对成功，写入完整信息
                    self.reporter.add_success(
                        group_id, filename1, filename2, entry1, entry2,
                        r1, r2, nmatch, rmsd_or_reason
                    )

            except Exception as e:
                # refcode无效（尝试提取已成功读取的entry的晶胞参数）
                self.reporter.add_failure(
                    group_id, filename1, filename2, r1, r2,
                    reason=f"refcode处理失败: {str(e)}",
                    entry1=entry1, entry2=entry2  # 可能有一个entry有效
                )

        print()  # 换行避免进度覆盖


def main():
    input_csv = "3.smiles_identifier_matched_preprocessed.csv"
    detailed_file = "detailed_pair_results_with_cell_params.csv"  # 含晶胞参数的详细文件
    summary_file = "filename_isomorphism_summary.csv"  # 新增identifier1和identifier2的汇总文件

    print("正在读取输入CSV并整理数据...")
    try:
        encodings = ['gbk', 'utf-8', 'latin-1', 'utf-16']
        for encoding in encodings:
            try:
                df = pd.read_csv(input_csv, encoding=encoding)
                print(f"✅ 成功使用编码 {encoding} 读取CSV")
                break
            except UnicodeDecodeError:
                continue
        else:
            raise Exception("❌ 所有编码尝试失败")
    except Exception as e:
        print(f"❌ 读取CSV错误: {e}")
        return

    # 列索引定义（关键：新增第四列和第八列的identifier）
    group_id_col = 0  # 第1列：group_id
    filename1_col = 1  # 第2列：filename1
    identifier1_col = 3  # 第4列：identifier1（需新增到汇总文件）
    filename2_col = 5  # 第6列：filename2
    identifier2_col = 7  # 第8列：identifier2（需新增到汇总文件）
    id1_col = 3  # 第4列：identifier1（refcode列表，用于生成比对对）
    id2_col = 7  # 第8列：identifier2（refcode列表，用于生成比对对）

    # 整理数据：
    # pairs_with_meta：用于比对的refcode对
    # all_group_info：用于汇总的行级信息（包含identifier1和identifier2）
    pairs_with_meta = []
    all_group_info = {}  # 格式：group_id -> (filename1, identifier1, filename2, identifier2)

    for _, row in df.iterrows():
        group_id = row.iloc[group_id_col]
        filename1 = row.iloc[filename1_col]
        identifier1 = row.iloc[identifier1_col]  # 提取第四列：identifier1
        filename2 = row.iloc[filename2_col]
        identifier2 = row.iloc[identifier2_col]  # 提取第八列：identifier2

        # 存储当前行的完整信息（用于汇总文件）
        all_group_info[group_id] = (filename1, identifier1, filename2, identifier2)

        # 分割refcode列表（用于生成比对对）
        r1_list = [r.strip() for r in str(row.iloc[id1_col]).split(',') if r.strip()]
        r2_list = [r.strip() for r in str(row.iloc[id2_col]).split(',') if r.strip()]

        # 生成所有refcode对
        for r1 in r1_list:
            for r2 in r2_list:
                pairs_with_meta.append((group_id, filename1, filename2, r1, r2))

    print(f"✅ 共生成 {len(pairs_with_meta)} 对refcode，覆盖 {len(all_group_info)} 行数据")

    # 开始比对
    print("\n开始比对（含晶胞参数实时写入）...")
    reporter = RealTimeReporter(detailed_file_path=detailed_file)
    analysis = PairAnalysis(pairs_with_meta=pairs_with_meta, reporter=reporter)

    try:
        analysis.run()
    except Exception as e:
        print(f"\n❌ 比对中断: {e}")
        return

    # 生成包含identifier1和identifier2的总结果
    print("\n✅ 比对完成，生成最终总结果...")
    summary = reporter.get_summary(all_group_info)
    with open(summary_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        # 表头新增“identifier1”和“identifier2”
        writer.writerow([
            "group_id",
            "filename1", "identifier1",  # filename1后紧跟identifier1（第四列）
            "filename2", "identifier2",  # filename2后紧跟identifier2（第八列）
            "同构状态（是/否）"
        ])
        for row in summary:
            writer.writerow(row)

    print("\n📁 结果文件：")
    print(f"   - 详细结果（含晶胞参数）：{detailed_file}")
    print(f"   - 最终总结果（含identifier）：{summary_file}")


if __name__ == "__main__":
    main()