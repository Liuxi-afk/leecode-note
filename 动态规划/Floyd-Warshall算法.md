解决**所有节点对之间最短路径（APSP）** 的经典算法
Floyd–Warshall 算法用的是一种 按允许中间顶点数量逐步扩展 的动态规划思想：定义 $dp[k][i][j]$ 为“只允许通过编号 ≤ k 的顶点作为中间顶点时，从 i 到 j 的最短路径长度”，并用递推式

$$dp[k][i][j] = \min(dp[k-1][i][j],\ dp[k-1][i][k] + dp[k-1][k][j])$$

逐步把可用中间顶点增加到所有顶点。时间复杂度 $O(n^3)$，空间 $O(n^2)$

## 三要素
### 状态含义
> 把“路径可以使用哪些中间点”作为状态维度。
- 状态：令 $dp[k][i][j]$ 表示:**从顶点i到顶点j的最短路径距离,且路径中间节点只能来自集合{1,2,...,k}.(起点i和终点j可以是任意编号,不受限制)**
- 目标:球 $dp[n][i][j]$(n为顶点数),即允许任意顶点作为中间点时的最短距离.

这个状态设计把“路径的复杂结构”转化为“中间点集合大小”的逐步扩张，方便用递推构造更复杂的解。

### 状态转移
考虑允许的中间点集合从{1,...,k - 1}扩展到{1,...,k},任意从i到j的最短路径要么:
1. 不经过顶点k作为中间点 --> 最短路径和 $dp[k][i][j]$ 相同
2. 经过顶点k(至少一次) --> 那么可以把路径拆成 i --> k 与k --> j,而这两段中间点都只允许来自{1,2,...,k - 1},其代价为 $dp[k - 1][i][k] + dp[k - 1][k][j]$

因此:

$$dp[k][i][j] = \min(dp[k-1][i][j],\ dp[k-1][i][k] + dp[k-1][k][j])$$

### 初始条件(边界)
- $dp[0][i][j] = $
  - 0 当 $i = j$;
  - $w[i][j]$如果存在边i --> j且权重为 $w[i][j]$
  - $\infty$ 如果 $i \neq j $且没有直接边
相当于初始化邻接矩阵(不存在的边设为 $\infty$)

## 算法框架
### 构建dp矩阵
```python
for i in range(n):
  for j in range(n):
    dp[i][j] = weight[i][j] #先初始化weight矩阵,对图中边记录临近权重,无边则记为float('inf')
for k in range(n):
  for i in range(n):
    for j in range(n):
      dp[i][j] = min(dp[i][j],dp[i][k] + dp[k][j])
```
### 路径重建（如何找出具体路径）
 常用做法: 维护一个 $next[i][j]$ 或 $pre[i][j]$ 矩阵.  
 $next[i][j]$表示表示从 i 到 j 的最短路径上的第一步应该走到哪个顶点。有了 $next$,从i到j的完整路径可以通过反复查询 $next$得到.
- 初始化:
  - 如果存在边i --> j ,令 $next[i][j] = j$ (或 $pre[j][i] = i$);否则 $None$.
  - $next[i][i] = i$(可选)
- 在内层更新时,如果 $dp[i][k] + dp[k][j] < dp[i][j]$,同时更新$next[i][j] = next[i][k]$ (表示要从i到j,第一步应走向到i -- > k路径的第一跳).
- 处理负权与负环
  - Floyd–Warshall 支持 负边权（只要没有负权环），因为递推式没有前提“非负”。
  - 若存在某个i使得最后 $dp[i][i] < 0$ ，说明存在一个能从i回到i且总权重为负的负环（可达到并绕行）。
    - 一旦有负环，最短路径在意义上就不存在（可以无限减小）。在算法中检测到 $dp[i][i] < 0$可判定存在负环。
  - 若图中存在负环，路径重建要小心：有负环的路径可能被“无限改进”，标准的 Floyd–Warshall 在检测到负环后通常会标注哪些对对之间的最短路径不再有界。
> 初始化dp矩阵,更新dp矩阵和next矩阵(从i到j最短路径的第一步选择),以及检测负环
```python
for i in range(n):
  for j in range(n):
    dp[i][j] = weight[i][j] #先初始化weight矩阵,对图中边记录临近权重,无边则记为float('inf')
for k in range(n):
  for i in range(n):
    for j in range(n):
      newdist = dp[i][k] + dp[k][j]
      if newdist < dp[i][j]:
        dp[i][j] = newdist
        next[i][j] = next[i][k] #更新最短路径第一步
# 检测负环（任何 dist[i][i] < 0 表示负环）
neg_cycle_nodes = [i for i in range(n) if dist[i][i] < 0]
```
- 重建路径 $path[i][j]$ :若 $next[i][j]$ 为 $None$,则无路径;否则从 $u = i$开始,重复 $u = next[u][j]$直到 $u == j$ ,收集顶点即为路径.
```python
def reconstruct_path(nxt, u, v):
    """
    使用 nxt 矩阵重建从 u 到 v 的路径（顶点序列）。
    如果不存在路径，返回 [] 或 None。
    注意：在使用前应先检查是否存在负环影响（dist 中 dist[x][x] < 0）。
    """
    if nxt[u][v] is None:
        return []  # 无路径

    path = [u]
    cur = u
    # 防止异常无限循环：若 nxt[cur][v] 为 None 则无路径
    while cur != v:
        cur = nxt[cur][v]
        if cur is None:
            return []  # 不可达（保险性检查）
        path.append(cur)
        # 若路径异常过长可额外检测循环（例如 len(path) > n），但若有负环需提前检测
    return path
```
