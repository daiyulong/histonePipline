#!/usr/bin/python3
# -*- coding:utf-8 -*-
# @Time   : 2022/11/3
# @Author : Dai Yulong

import sys
import re
import os
import argparse
import platform

"""Data Format Conversion"""
def run_data_conversion(inFile, pn, bin, info,olog):
    if not os.path.exists(pn):
        os.mkdir(pn)
    if os.path.exists(inFile):
        # 格式转换
        print("### Data Format Conversio ###")
        olog.write(
            "Rscript {0}/preAnalysisHistone.R {1} {2}\n".format(bin, inFile, pn))
        os.system(
            "Rscript {0}/preAnalysisHistone.R {1} {2}\n".format(bin, inFile, pn))
        # 参数读取，生成group.txt 和 input.txt文件
        olog.write("Rscript {}/readParaHistone.R {} {} {}\n".format(bin, info, pn + "/infile.csv", pn))
        os.system("Rscript {}/readParaHistone.R {} {} {}".format(bin, info, pn + "/infile.csv", pn))


"""Sample Analysis"""
def run_sample_analysis(groupFile, scaleMethod, pn, bin, olog, delmin):
    of = pn + "/SampleAnalysis"
    if not os.path.exists(of):
        os.mkdir(of)
    if os.path.exists(groupFile):
        # 探索性绘图
        print("\n### Run Sample Analysis ###\n")
        olog.write(
            "Rscript {0}/matrixDrawHistone.R {1} {2} {3} Intensity {4} {5}\n".format(
                bin, pn + "/infile.csv", groupFile, of, scaleMethod, delmin))
        os.system(
            "Rscript {0}/matrixDrawHistone.R {1} {2} {3} Intensity {4} {5}".format(
                bin, pn + "/infile.csv", groupFile, of, scaleMethod, delmin))

        olog.write("Rscript {}/svdiPCAHistone.R {} {} {} uv\n".format(bin, pn + "/infile.csv", groupFile, of))
        os.system("Rscript {}/svdiPCAHistone.R {} {} {} uv".format(bin, pn + "/infile.csv", groupFile, of))



"""统计检验差异筛选"""
def run_stat_analysis(L, bin, pn, olog, fc=1.5, pvalue=0.05):
    print("\n### Run Statistic Analysis ###\n")
    for comp in L:
        of = pn + '/' + comp
        infile = pn + '/' + comp + "/input.txt"
        groupfile = pn + '/' + comp + '/group.txt'
        if os.path.exists(infile) and os.path.exists(groupfile):
            
            print("\n## Run " + comp + " ##\n")
            # 差异筛选
            if L[comp] == 'ttest':
                olog.write("Rscript {}/multitest.R {} {} {} stt {} {}\n".format(bin, infile, groupfile, pn + '/' + comp, fc, pvalue))
                os.system("Rscript {}/multitest.R {} {} {} stt {} {}".format(bin, infile, groupfile, pn + '/' + comp, fc, pvalue))
            elif L[comp] == 'pairedttest':
                olog.write("Rscript {}/multitest.R {} {} {} ptt {} {} TRUE\n".format(bin, infile, groupfile, pn + '/' + comp, fc, pvalue))
                os.system("Rscript {}/multitest.R {} {} {} ptt {} {} TRUE".format(bin, infile, groupfile, pn + '/' + comp, fc, pvalue))

            #绘制堆叠图
            olog.write("Rscript {}/drawStackedBar.R {} {} {} {}\n".format(bin, infile, pn + '/infile.csv', groupfile, pn + '/' + comp))
            os.system("Rscript {}/drawStackedBar.R {} {} {} {}".format(bin, infile, pn + '/infile.csv', groupfile, pn + '/' + comp))

            #绘制errorBar图
            olog.write("Rscript {}/drawerrorBar.R {} {}/test.txt {} {}\n".format(bin, pn + '/infile.csv', of, groupfile, of))
            os.system("Rscript {}/drawerrorBar.R {} {}/test.txt {} {}".format(bin,  pn + '/infile.csv', of, groupfile, of))

            #绘制分组heatmap图
            olog.write("Rscript {}/drawHeatmapHistone.R {} {} {}\n".format(bin, infile, groupfile, of))
            os.system("Rscript {}/drawHeatmapHistone.R {} {} {}".format(bin,  infile, groupfile, of))




"""获得比较组信息"""
def get_comparison(infile):
    ifile = open(infile)
    L = dict()
    for line in ifile.readlines()[1:]:
        row = line.rstrip().split('\t')
        row[0] = row[0].lower()
        if row[0] == 'ttest' or row[0] == 'pairedttest':
            L["{}.vs.{}".format(row[1], row[2])] = row[0]
    return (L)


"""创建一个目录"""
def createdir(d):
    if not os.path.exists(d):
        os.makedirs(d)

"""获取当前程序所在的目录"""
def get_absdir():
    if not sys.argv[0].endswith('.py') and not sys.argv[0].endswith('.pyw'):
        absDir = os.path.dirname(sys.executable)
    else:
        absDir = os.path.dirname(os.path.abspath(__file__))

    return absDir


def get_bin():
    if platform.system() == 'Windows':
        bin = get_absdir()
    else:
        bin = '/usr/local/bio/PAA/bin'
    return bin

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="HistoneWflow")
    parser.add_argument('-o', action='store', dest='pn', default='', help='output path')
    parser.add_argument('-hf', action='store', dest='hf', default='inputfile.txt', help='histone raw inputfile.txt ')
    parser.add_argument('-info', action='store', dest='info', default='parameter.xlsx', help='input parameters.xlsx')
    parser.add_argument('-fc', action='store', dest='fc', default=1.5, help='fold change cutoff, default is 1.5')
    parser.add_argument('-pvalue', action='store', dest='pvalue', default=0.05, help='p-value cutoff, default is 0.05')
    parser.add_argument('-sm', action='store', dest='scaleMethod', default='log2', help='scale method for data draw')
    parser.add_argument('-mm', action='store', dest='maskminvalue', default=False,
                        help="是否屏蔽掉最小值后再做箱线图，默认是False")

    p = parser.parse_args()

    if not os.path.exists(p.hf) or not os.path.exists(p.info) or p.pn == '':
        parser.print_help()
        sys.exit()

    bin = get_bin()
    createdir(p.pn)
    olog = open(p.pn + '/command.log', 'a')

    ## 表格格式转换
    run_data_conversion(p.hf, p.pn, bin, p.info, olog)


    ## 获得比较组信息
    L = get_comparison(p.pn + '/comparison.txt')
    print(L)

    ## SampleAnalysis
    del_min = ''
    if p.maskminvalue is False:
        del_min = 'FALSE'
    else:
        del_min = 'TRUE'
    run_sample_analysis(p.pn+'/groupFile.txt', p.scaleMethod, p.pn, bin, olog, del_min)

    ## 统计分析
    run_stat_analysis(L, bin, p.pn, olog, p.fc, p.pvalue)


    print("完成分析。")
    olog.close()
