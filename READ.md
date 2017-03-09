# 说明

`/src` 中为项目代码，其中`hyperloglog.cc`和`LFM.cc`分别是两个代码的源文件

项目用到了boost库的部分内容，[下载地址](http://www.boost.org)。

`/test script` 下是测试脚本，用到了 Python 的`numpy`和`pandas`库，生成的数据条目满足zipf分布。

你也可以直接使用：


```
python zipf-gen-card.py 100000
```

表示使用脚本文件生成100000条数据并将其保存在文件中。



```
ELDA.exe 0
```

表示使用LFM进行势分析

```
ELDA.exe 1
```

表示使用Hyperloglog进行势分析
