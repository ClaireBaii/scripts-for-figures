

## A. 明天开工前 5 分钟：确认终端与代理状态

1. 打开 **Antigravity → Terminal**，确认默认终端是 **PowerShell**（建议）
2. 先测网络是否能直连 GitHub（有 VPN 也可能不需要额外代理配置）：

```powershell
ping github.com
```

如果 ping 被禁不代表不能用 git；关键用下面这个测试：

```powershell
git ls-remote https://github.com/ClaireBaii/scripts-for-figures.git
```

> 这条命令成功返回一堆 hash/refs，就说明：网络 + 认证（至少读权限）没问题。

---

## B. 安装 Git（让 Antigravity 的 Source Control 能提交/拉取/推送）

### B1) 安装 Git for Windows

如果你还没装 Git：

```powershell
winget search git
winget install Git.Git
```

装完重开 Antigravity（让它重新检测 git），然后检查：

```powershell
git --version
where git
```

### B2) Git 全局配置（只做一次）

```powershell
git config --global user.name  "Jacob"
git config --global user.email "你的邮箱"

git config --global core.autocrlf input
git config --global core.longpaths true
```

### B3)（可选）若 git 访问 GitHub 卡住：给 git 单独走 Clash 端口

Clash Verge 常见本地端口是 `127.0.0.1:7890`（以你本机为准）。只有在 `git ls-remote` 失败时才做：

```powershell
git config --global http.proxy  http://127.0.0.1:7890
git config --global https.proxy http://127.0.0.1:7890
```

恢复（不用代理时）：

```powershell
git config --global --unset http.proxy
git config --global --unset https.proxy
```

---

## C. Clone 仓库到本地 + 让 Antigravity 关联仓库

### C1) 选择工作目录并拉取

建议放到路径短的位置，例如 `D:\work`：

```powershell
mkdir D:\work -Force
cd D:\work

git clone https://github.com/ClaireBaii/scripts-for-figures.git
cd scripts-for-figures
```

### C2) 用 Antigravity 打开并验证 Git 集成

* Antigravity → **Open Folder** → 选择 `D:\work\scripts-for-figures`
* 左侧 **Source Control** 应该能看到：

  * 当前分支
  * 变更列表
  * commit / sync 按钮

如果 Source Control 提示找不到 Git：
到设置里把 git 路径指向（通常）：
`C:\Program Files\Git\cmd\git.exe`

---

## D. 安装 Conda（推荐 Miniforge）并用 environment.yml 配 R 环境

### D1) 安装 Miniforge（推荐）或 Miniconda

用 winget 搜索（不同机器包名可能不同）：

```powershell
winget search miniforge
winget search miniconda
```

如果你找不到 Miniforge，就装 Miniconda（也可以）。装完后检查：

```powershell
conda --version
```

### D2) 初始化 conda 到 PowerShell（关键）

```powershell
conda init powershell
```

然后 **彻底关闭 Antigravity / 终端再重开**。

如果报执行策略问题（常见）：

```powershell
Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned
```

---

## E. 在项目里创建 conda 环境并确认 R 可用

### E1) 找到 environment.yml 并创建环境

在仓库根目录：

```powershell
cd D:\work\scripts-for-figures
dir
```

如果看到 `environment.yml`：

```powershell
conda env create -f environment.yml
conda env list
```

激活环境（环境名以 yml 的 `name:` 为准）：

```powershell
conda activate <env_name>
```

验证 R：

```powershell
R --version
where R
Rscript --version
```

### E2)（可选）环境已存在就更新（更稳）

```powershell
conda env update -f environment.yml --prune
```

### E3)（可选）若 conda 下载慢/失败：让 conda 走 Clash

只在 conda 明显连不上时设置：

```powershell
conda config --set proxy_servers.http  http://127.0.0.1:7890
conda config --set proxy_servers.https http://127.0.0.1:7890
```

取消：

```powershell
conda config --remove-key proxy_servers
```

---

## F. 让 Antigravity “用 conda 环境里的 R” 跑脚本

### F1) 安装 R 扩展

在 Antigravity 的 Extensions 里装：

* **R**（语言支持）

### F2) 配置 R.exe 路径（最关键的一步）

在**已激活 conda 环境**的终端里执行：

```powershell
where R
```

你会得到类似：
`C:\Miniforge3\envs\<env_name>\Library\bin\R.exe`

把这个路径写入 Antigravity Settings（JSON）里（通常与 VS Code 兼容）：

```json
{
  "r.rterm.windows": "C:\\Miniforge3\\envs\\<env_name>\\Library\\bin\\R.exe"
}
```

> 如果 Antigravity 的 R 扩展设置项名称与 VS Code 不完全一致：就用设置搜索关键词 `rterm` / `R.exe`，把路径填进去。核心就是“指向 conda env 的 R.exe”。

### F3) 推荐运行方式（最少踩坑）

在 Antigravity 终端里：

```powershell
conda activate <env_name>
cd D:\work\scripts-for-figures
Rscript path\to\your_script.R
```

---

## G. 明天你要“改图代码优化效果”：建议按 Git 分支工作，随时可回滚

### G1) 开一个专用分支（强烈建议）

```powershell
cd D:\work\scripts-for-figures
git checkout -b improve-fig-style
```

### G2) 建立“基线输出”（改代码前先跑一遍）

* 先跑一遍脚本，保存输出图片到一个明确目录（例如 `outputs_baseline/`，目录名以项目实际为准）
* 之后每次 AI 改代码 → 重新跑 → 对比 baseline

### G3) 每次改动小步提交

```powershell
git status
git add -A
git commit -m "Improve figure aesthetics: <简述>"
```

---

## H. 最小自检：确认“R 能在 conda 环境里生成图片”

如果你想不依赖项目脚本，先做一次最小测试（排除图形设备问题）：

```powershell
conda activate <env_name>
Rscript -e "png('test.png', width=2000, height=1500, res=300); plot(1:10); dev.off()"
dir test.png
```

看到 `test.png` 就说明：R 环境 + 图片输出链路没问题。

---

## 你明天的“流水线”按这个顺序执行就行

1. 安装 Git（若未装）→ 配 git config
2. `git ls-remote` 测仓库可达
3. `git clone` 到 `D:\work`
4. 安装 conda（Miniforge/Miniconda）→ `conda init powershell` → 重开
5. `conda env create -f environment.yml` → `conda activate` → `where R`
6. Antigravity 配 `r.rterm.windows` 指向 conda env 的 `R.exe`
7. 终端里 `Rscript xxx.R` 跑项目出图
8. 开分支 `improve-fig-style` → 用 AI 小步改 → 跑 → 对比 → commit

