# Differential-evolution
Differential evolution

這是一份群體智慧的期末報告：
使用差分進化演算法 (Differential evolution, DE)與其三個DE變型體(APTDE, EPSDE, MPEDE)去找十個目標函數 (objective function) 的全域最優解 (global optimal solution)。APTDE採用自適應群體數；EPSED採用自適應參數；MPEDE採用多群+自適應參數。整體表現上，（優）MPEDE > EPSDE > APTDE > Basic DE（差）。

報告全文：
1. 差分進化演算法 Differential evolution - 1
2. https://tsungsquare.blogspot.com/2021/06/differential-evolution-1.html
3. 差分進化演算法 Differential evolution - 2｜APTDE｜自適應群體
4. https://tsungsquare.blogspot.com/2021/06/differential-evolution-2aptde.html
5. 差分進化演算法 Differential evolution - 3｜EPSDE、MPEDE｜自適應參數、多群
6. https://tsungsquare.blogspot.com/2021/06/differential-evolution-3epsdempede.html
7. 差分進化演算法 Differential evolution - 4｜APTDE, EPSDE, MPEDE比較結果
8. https://tsungsquare.blogspot.com/2021/06/differential-evolution-4aptde-epsde.html

參考文獻：
DE
Storn, R., & Price, K. (1997). Differential evolution - A simple and efficient heuristic for global optimization over continuous spaces. Journal of Global Optimization, 11(4), 341-359. https://doi.org/Doi/10.1023/A:1008202821328

EPSDE
Mallipeddi, R., Suganthan, P. N., Pan, Q. K., & Tasgetiren, M. F. (2011). Differential evolution algorithm with ensemble of parameters and mutation strategies. Applied Soft Computing, 11(2), 1679-1696. https://doi.org/10.1016/j.asoc.2010.04.024

APTDE
Zhu, W., Tang, Y., Fang, J.-a., & Zhang, W. (2013). Adaptive population tuning scheme for differential evolution. Information Sciences, 223, 164-191. https://doi.org/10.1016/j.ins.2012.09.019 

MPEDE
Wu, G., Mallipeddi, R., Suganthan, P. N., Wang, R., & Chen, H. (2016). Differential evolution with multi-population based ensemble of mutation strategies. Information Sciences, 329, 329-345. https://doi.org/10.1016/j.ins.2015.09.009
