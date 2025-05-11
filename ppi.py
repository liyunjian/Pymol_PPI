import pymol
from pymol import cmd
import pandas as pd
import numpy as np
import os

def calculate_interactions(target_selection="sele", output_file="interactions.xlsx"):
    """
    计算选定分子与其他分子之间的相互作用，并将结果导出到Excel表格
    
    参数:
        target_selection: 目标分子的选择表达式
        output_file: 输出的Excel文件名
    """
    # 确保PyMOL已初始化
    if not hasattr(pymol, 'finish_launching'):
        pymol.finish_launching()
    
    # 获取目标选择的原子
    target_atoms = cmd.get_model(target_selection)
    if len(target_atoms.atom) == 0:
        print("错误：未找到目标分子，请确保选择了正确的分子")
        return
    
    # 获取所有其他原子
    all_atoms = cmd.get_model("all and not " + target_selection)
    
    # 准备存储相互作用的列表
    interactions = []
    
    # 遍历目标分子的每个原子
    for target_atom in target_atoms.atom:
        target_coord = np.array([target_atom.coord[0], target_atom.coord[1], target_atom.coord[2]])
        target_resn = target_atom.resn  # 残基名称
        target_resi = target_atom.resi  # 残基编号
        target_chain = target_atom.chain  # 链ID
        target_name = target_atom.name  # 原子名称
        target_elem = get_element(target_atom)  # 元素
        
        # 遍历所有其他原子
        for other_atom in all_atoms.atom:
            other_coord = np.array([other_atom.coord[0], other_atom.coord[1], other_atom.coord[2]])
            other_elem = get_element(other_atom)  # 元素
            
            # 计算距离
            distance = np.linalg.norm(target_coord - other_coord)
            
            # 根据不同相互作用类型的距离阈值检测相互作用
            interaction_type = determine_interaction_type(target_atom, other_atom, distance)
            
            # 如果检测到相互作用，则添加到列表
            if interaction_type:
                interactions.append({
                    "目标残基": f"{target_resn}{target_resi}",
                    "目标链": target_chain,
                    "目标原子": target_name,
                    "目标元素": target_elem,
                    "相互作用残基": f"{other_atom.resn}{other_atom.resi}",
                    "相互作用链": other_atom.chain,
                    "相互作用原子": other_atom.name,
                    "相互作用元素": other_elem,
                    "距离(Å)": round(distance, 2),
                    "相互作用类型": interaction_type
                })
    
    # 如果没有找到相互作用
    if not interactions:
        print("未找到相互作用")
        return
    
    # 创建DataFrame并导出到Excel
    df = pd.DataFrame(interactions)
    
    # 按照目标残基、相互作用残基和距离排序
    df = df.sort_values(by=["目标链", "目标残基", "相互作用链", "相互作用残基", "距离(Å)"])
    
    # 导出到Excel
    df.to_excel(output_file, index=False)
    print(f"相互作用已导出到 {os.path.abspath(output_file)}")

def determine_interaction_type(atom1, atom2, distance):
    """
    根据原子类型和距离确定相互作用类型
    
    根据不同相互作用类型的特定距离阈值判断：
    - 氢键: 3.5 Å以内
    - 疏水相互作用: 5.0 Å以内
    - 盐桥: 4.0 Å以内
    - π-π相互作用: 6.0 Å以内
    - 离子相互作用: 6.0 Å以内
    - 水桥: 3.5 Å以内
    - 范德华力: 4.0 Å以内
    """
    # 获取原子元素
    elem1 = get_element(atom1)
    elem2 = get_element(atom2)
    
    # 氢键检测 (N-O, O-N, O-O 等) - 距离3.5 Å以内
    if ((elem1 in ['N', 'O'] and elem2 in ['N', 'O']) and distance <= 3.5):
        return "氢键"
    
    # 盐桥检测 (带电荷的原子对) - 距离4.0 Å以内
    if is_charged(atom1) and is_charged(atom2) and distance <= 4.0:
        return "盐桥"
    
    # π-π堆积 (芳香环之间) - 距离6.0 Å以内
    if is_aromatic(atom1) and is_aromatic(atom2) and distance <= 6.0:
        return "π-π相互作用"
    
    # 疏水相互作用 - 距离5.0 Å以内
    if is_hydrophobic(elem1) and is_hydrophobic(elem2) and distance <= 5.0:
        return "疏水相互作用"
    
    # 离子相互作用 - 距离6.0 Å以内
    if is_ionic(atom1) and is_ionic(atom2) and distance <= 6.0:
        return "离子相互作用"
    
    # 水桥检测 - 距离3.5 Å以内
    if is_water(atom1) and (elem2 in ['N', 'O']) and distance <= 3.5:
        return "水桥"
    if is_water(atom2) and (elem1 in ['N', 'O']) and distance <= 3.5:
        return "水桥"
    
    # 范德华力相互作用 - 距离4.0 Å以内
    if distance <= 4.0:
        return "范德华力"
    
    # 如果距离大于所有阈值，则不返回相互作用类型
    return None

def get_element(atom):
    """
    从原子名称中提取元素符号
    """
    # 尝试使用elem属性（某些PyMOL版本支持）
    if hasattr(atom, 'elem') and atom.elem.strip():
        return atom.elem.strip()
    
    # 否则从原子名称中提取元素符号（通常是名称的第一个或前两个字符）
    name = atom.name.strip()
    
    # 处理常见的原子命名约定
    if name.startswith(('C', 'N', 'O', 'S', 'P', 'H', 'F', 'I', 'K')):
        return name[0]
    elif name.startswith(('CL', 'BR', 'FE', 'MG', 'ZN', 'CA')):
        return name[:2].capitalize()
    else:
        # 默认返回第一个字符
        return name[0] if name else 'X'

def is_charged(atom):
    """检查原子是否带电荷"""
    charged_atoms = {
        'ARG': ['NH1', 'NH2', 'NE'],
        'LYS': ['NZ'],
        'ASP': ['OD1', 'OD2'],
        'GLU': ['OE1', 'OE2'],
        'HIS': ['ND1', 'NE2']
    }
    
    if atom.resn in charged_atoms and atom.name in charged_atoms[atom.resn]:
        return True
    return False

def is_ionic(atom):
    """检查原子是否参与离子相互作用"""
    # 离子相互作用通常涉及带电荷的原子，与is_charged类似但可能有所不同
    ionic_residues = ['ARG', 'LYS', 'ASP', 'GLU', 'HIS']
    ionic_atoms = {
        'ARG': ['NH1', 'NH2', 'NE'],
        'LYS': ['NZ'],
        'ASP': ['OD1', 'OD2'],
        'GLU': ['OE1', 'OE2'],
        'HIS': ['ND1', 'NE2']
    }
    
    if atom.resn in ionic_residues and atom.name in ionic_atoms.get(atom.resn, []):
        return True
    return False

def is_aromatic(atom):
    """检查原子是否属于芳香环"""
    aromatic_residues = ['PHE', 'TYR', 'TRP', 'HIS']
    aromatic_atoms = {
        'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'TRP': ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
        'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2']
    }
    
    if atom.resn in aromatic_residues and atom.name in aromatic_atoms.get(atom.resn, []):
        return True
    return False

def is_hydrophobic(element):
    """检查元素是否疏水"""
    # 碳原子通常被认为是疏水的
    hydrophobic_elements = ['C']
    return element in hydrophobic_elements

def is_water(atom):
    """检查原子是否属于水分子"""
    # 水分子通常被标记为HOH, WAT, H2O等
    water_residues = ['HOH', 'WAT', 'H2O', 'TIP', 'TIP3', 'TIP4']
    return atom.resn in water_residues

def load_and_analyze(pdb_file, selection_expression, output_file="interactions.xlsx"):
    """
    加载PDB文件并分析相互作用
    
    参数:
        pdb_file: PDB文件路径
        selection_expression: 选择表达式，例如 "chain A"
        output_file: 输出文件名
    """
    # 加载PDB文件
    cmd.load(pdb_file)
    
    # 选择目标分子
    cmd.select("target_selection", selection_expression)
    
    # 计算相互作用
    calculate_interactions("target_selection", output_file)

def analyze_chains_interaction(pdb_file, target_chains, output_file="chain_interactions.xlsx"):
    """
    分析指定链组与其他链之间的相互作用
    
    参数:
        pdb_file: PDB文件路径
        target_chains: 目标链的列表，例如 ["A", "B"]
        output_file: 输出文件名
    """
    # 加载PDB文件
    cmd.load(pdb_file)
    
    # 构建选择表达式
    chain_selection = " or ".join([f"chain {chain}" for chain in target_chains])
    
    # 选择目标链
    cmd.select("target_chains", chain_selection)
    
    # 计算相互作用
    calculate_interactions("target_chains", output_file)
    
    print(f"已分析链 {', '.join(target_chains)} 与其他链之间的相互作用")

# 将函数扩展到PyMOL命令行
cmd.extend("calculate_interactions", calculate_interactions)
cmd.extend("load_and_analyze", load_and_analyze)
cmd.extend("analyze_chains_interaction", analyze_chains_interaction)

# 使用示例
print("PyMOL 分子相互作用分析脚本已加载")
print("使用方法:")
print("1. 在PyMOL中选择目标分子")
print("2. 运行: calculate_interactions selection_name, output_file.xlsx")
print("   例如: calculate_interactions chain A, interactions.xlsx")
print("或者:")
print("3. 直接加载PDB并分析: load_and_analyze pdb_file, selection_expression, output_file.xlsx")
print("   例如: load_and_analyze 1abc.pdb, chain A, interactions.xlsx")
print("或者:")
print("4. 分析多条链与其他链的相互作用: analyze_chains_interaction pdb_file, [\"A\", \"B\"], output_file.xlsx")
print("   例如: analyze_chains_interaction 1abc.pdb, [\"A\", \"B\"], AB_interactions.xlsx")
print("\n相互作用类型及其距离阈值:")
print("- 氢键: 3.5 Å以内")
print("- 疏水相互作用: 5.0 Å以内")
print("- 盐桥: 4.0 Å以内")
print("- π-π相互作用: 6.0 Å以内")
print("- 离子相互作用: 6.0 Å以内")
print("- 水桥: 3.5 Å以内")
print("- 范德华力: 4.0 Å以内")

"""
- 启动PyMOL
- 运行脚本：在PyMOL命令行中输入 run PPI.py
- 加载您的PDB文件： load your_structure.pdb
- 选择目标分子：例如 select target, chain A+B
- 分析相互作用： calculate_interactions target, interactions.xlsx
"""