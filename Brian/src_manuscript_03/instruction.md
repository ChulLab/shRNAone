## part A
1. 清洗湿实验结果
2. 使用预测软件对湿实验数据进行预测，得到预测结果。并统一坐标加以清洗
3. 处理一下offset
4. 与湿实验结果进行相关性分析，发现casRX表现尚可，但是无法合理解释shRNAone 和 psp
## part B
1. 使用1nt尺度，遍历W,L,U组合，基于TTR扫描全局。设置criteria提取合理显著组合。观察cor分布状况，确定使用3prime
2. 固定WLU的两个参数，改变另一个参数。得到U是影响cor最大的参数
## part C
1. 固定W=80，L=40，去遍历U，得到不同工具下最佳的U和合理的U范围。
2. 从W=80，L=40中得到各个工具中U的稳定平台，限定这些U，然后寻找一组W,L能够在平台上最稳定的去解释三个工具。此处得到基于U平台的最佳参数（组合）
3. 前两者均能顺利通过siRNA ground truth的验证
## part D
1. 确定好基于W=80,L=40的一组最佳参数组合，以及基于WL screen得到的最佳参数组合 (W35L20)
2. 两者均能顺利通过siRNA ground truth的验证 （shRNAone）
3. 两者均能/部分能解释PCSK9的湿实验结果。（shRNAone）
3. 能够良好区别TTR 的permutation定义的区域 or top percent定义的区域，对PCSK9也许也有一定的区别作用
4. 将两个参数与预测工具来源数据进行大横评
## part E
1. 使用permutation方法定义shRNAone-TTR, PCSK9优良区域和非优良区域。也许可以遍历不同rolling window
2. 观察区域之间accessibility分布
## part F
1. 使用这两个参数重新计算人类基因组数据
2. 为PSP提供指导序列