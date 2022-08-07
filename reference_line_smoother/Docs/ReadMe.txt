# optimization form
# min. 0.5*x^(T)*H*x + f^(T)*x
# s.t. l << A*x << u

# for osqp optimization, following website can be useful
# https://chowdera.com/2022/04/202204200603266904.html
# https://zhuanlan.zhihu.com/p/376122376
# https://robotology.github.io/osqp-eigen/md_pages_mpc.html

# for FemReferenceLineSmoother Algorithm, in following website has given related detail
# https://zhuanlan.zhihu.com/p/371585754
# https://zhuanlan.zhihu.com/p/523913529

# for SparseMatrix you can see here:
# https://blog.csdn.net/xuezhisdc/article/details/54633274
# https://eigen.tuxfamily.org/dox/group__QuickRefPage.html
# https://blog.csdn.net/qq_27806947/article/details/105545360

# here we should notice that:
## 1. Eigen::Solver can't be written as std::shared_ptr form, or it will make core dump
## 2. for FemReferenceLineSmoother, it input waypoints shouldn't be less than 3 points, 
##   or it win fail 


# log need to be change as Boost log