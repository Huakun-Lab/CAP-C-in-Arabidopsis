# Modified from https://github.com/ouyang-lab/CAPC/tree/master/src/filter_N.py
#!/python
# -- coding: utf-8 --
import gzip
import argparse


def parser_merged_reads_deal(Out_Ok1, Out_Ok2, Out_Discard1, Out_Discard2, t1, t2):
    for no, line in enumerate(zip(t1, t2)):
        if (no % 4) == 0:  # 第一行是fastq的
            R1_name = str(line[0], 'utf-8')  # 取第一个和第三个文件的序列名那一行
            R2_name = R1_name.replace('/1', '/2')

        elif (no % 4) == 1:  # 序列行（fastq的第二行）
            R1_seq = str(line[0], 'utf-8')
            R2_seq = str(line[1], 'utf-8')

            R1_length = len(line[0])
            R2_length = len(line[1])

            if R1_seq == R2_seq:
                file_handle = (Out_Discard1, Out_Discard2)
                print("%s\n%s\n+" % (R1_name, R1_seq), file=file_handle[0])
                print("%s\n%s\n+" % (R2_name, R2_seq), file=file_handle[1])
            else:

                if R1_length >= min_length and R2_length >= min_length:
                    file_handle = (Out_Ok1, Out_Ok2)
                    print("%s\n%s\n+" % (R1_name, R1_seq), file=file_handle[0])
                    print("%s\n%s\n+" % (R2_name, R2_seq), file=file_handle[1])
                else:
                    file_handle = (Out_Discard1, Out_Discard2)
                    print("%s\n%s\n+" % (R1_name, R1_seq), file=file_handle[0])
                    print("%s\n%s\n+" % (R2_name, R2_seq), file=file_handle[1])

        elif (no % 4) == 3:  # 质量行的处理
            quality_R1 = str(line[0], 'utf-8')
            quality_R2 = str(line[1], 'utf-8')

            print("%s" % quality_R1, file=file_handle[0])
            print("%s" % quality_R2, file=file_handle[1])


def parser_no_merge_reads_deal(Out_Ok1, Out_Ok2, Out_Short1, Out_Short2, Out_Discard1, Out_Discard2, In_y1, In_y2,
                               In_y3, In_y4):
    for no, line in enumerate(zip(In_y1, In_y2, In_y3, In_y4)):
        if (no % 4) == 0:  # 第一行是fastq的
            R1_name = str(line[0], 'utf-8')  # 取第一个和第三个文件的序列名那一行
            R2_name = str(line[2], 'utf-8')
        elif (no % 4) == 1:  # 序列行（fastq的第二行）

            select = None
            file_handle = None

            raw_seq_R1_P1_length = len(str(line[0], 'utf-8'))
            raw_seq_R1_P2_length = len(str(line[1], 'utf-8'))
            raw_seq_R2_P1_length = len(str(line[2], 'utf-8'))
            raw_seq_R2_P2_length = len(str(line[3], 'utf-8'))
            # 各序列的原始长度
            (l_R1_P1_length, l_R1_P2_length, l_R2_P1_length, l_R2_P2_length) = \
                (raw_seq_R1_P1_length - len(str(line[0], 'utf-8').lstrip("N")),
                 raw_seq_R1_P2_length - len(str(line[1], 'utf-8').lstrip("N")),
                 raw_seq_R2_P1_length - len(str(line[2], 'utf-8').lstrip("N")),
                 raw_seq_R2_P2_length - len(str(line[3], 'utf-8').lstrip("N")))
            # 记录序列开头处有几个N
            (r_R1_P1_length, r_R1_P2_length, r_R2_P1_length, r_R2_P2_length) = \
                (raw_seq_R1_P1_length - len(str(line[0], 'utf-8').rstrip("N")),
                 raw_seq_R1_P2_length - len(str(line[1], 'utf-8').rstrip("N")),
                 raw_seq_R2_P1_length - len(str(line[2], 'utf-8').rstrip("N")),
                 raw_seq_R2_P2_length - len(str(line[3], 'utf-8').rstrip("N")))
            # 记录序列末尾处有几个N

            seq_R1_P1_rmN = str(line[0], 'utf-8').strip("N")
            seq_R1_P2_rmN = str(line[1], 'utf-8').strip("N")
            seq_R2_P1_rmN = str(line[2], 'utf-8').strip("N")
            seq_R2_P2_rmN = str(line[3], 'utf-8').strip("N")

            four_seq_list = [str(line[0], 'utf-8'), str(line[1], 'utf-8'), str(line[2], 'utf-8'), str(line[3], 'utf-8')]
            # 将四个部分的序列数据存入数组，方便后续处理

            i = 0
            # 增加计数器，用以统计四个部分的序列数据中，哪些含有linker

            for seq_tmp in four_seq_list:
                if linker_seq in seq_tmp:
                    i += 1
            # 四个序列中，每有一个序列含有linker，则计数器加一，用以进行后续处理
            # 其中分为三种情况：
            # i == 4， i == 2, i == 0
            # 分别对应R1与R2中均含有linker，R1或R2中有一端含有linker，R1和R2中均不含有linker
            # 其中当i==4，这种情况中，R1和R2理论上为反向互补序列，但由于mismatch较多，在之前的seqprep处理时，无法将两个reads merge到一起
            # 所以出现了在R1和R2中均存在linker的现象。
            if i == 0:
                select = (1, 3)
                file_handle = (Out_Ok1, Out_Ok2)
                print("%s\n%s\n+" % (R1_name, seq_R1_P2_rmN), file=file_handle[0])
                print("%s\n%s\n+" % (R2_name, seq_R2_P2_rmN), file=file_handle[1])
            else:
                if len(seq_R1_P2_rmN) >= min_length and len(seq_R2_P2_rmN) >= min_length:
                    select = (1, 3)
                    file_handle = (Out_Ok1, Out_Ok2)
                    print("%s\n%s\n+" % (R1_name, seq_R1_P2_rmN), file=file_handle[0])
                    print("%s\n%s\n+" % (R2_name, seq_R2_P2_rmN), file=file_handle[1])
                elif len(seq_R1_P1_rmN) >= min_length and len(seq_R2_P2_rmN) >= min_length:
                    select = (0, 3)
                    file_handle = (Out_Short1, Out_Short2)
                    print("%s\n%s\n+" % (R1_name, seq_R1_P1_rmN), file=file_handle[0])
                    print("%s\n%s\n+" % (R2_name, seq_R2_P2_rmN), file=file_handle[1])
                elif len(seq_R1_P2_rmN) >= min_length and len(seq_R2_P1_rmN) >= min_length:
                    select = (1, 2)
                    file_handle = (Out_Short1, Out_Short2)
                    print("%s\n%s\n+" % (R1_name, seq_R1_P2_rmN), file=file_handle[0])
                    print("%s\n%s\n+" % (R2_name, seq_R2_P1_rmN), file=file_handle[1])
                elif len(seq_R1_P1_rmN) >= min_length and len(seq_R2_P1_rmN) >= min_length:
                    select = (0, 2)
                    file_handle = (Out_Short1, Out_Short2)
                    print("%s\n%s\n+" % (R1_name, seq_R1_P1_rmN), file=file_handle[0])
                    print("%s\n%s\n+" % (R2_name, seq_R2_P1_rmN), file=file_handle[1])
                else:
                    select = (1, 3)
                    file_handle = (Out_Discard1, Out_Discard2)
                    print("%s\n%s\n+" % (R1_name, seq_R1_P2_rmN), file=file_handle[0])
                    print("%s\n%s\n+" % (R2_name, seq_R2_P2_rmN), file=file_handle[1])

        elif (no % 4) == 3:  # 质量行的处理
            quality_R1_P1 = str(line[0], 'utf-8')[l_R1_P1_length:raw_seq_R1_P1_length - r_R1_P1_length]
            quality_R1_P2 = str(line[1], 'utf-8')[l_R1_P2_length:raw_seq_R1_P2_length - r_R1_P2_length]
            quality_R2_P1 = str(line[2], 'utf-8')[l_R2_P1_length:raw_seq_R2_P1_length - r_R2_P1_length]
            quality_R2_P2 = str(line[3], 'utf-8')[l_R2_P2_length:raw_seq_R2_P2_length - r_R2_P2_length]
            Quality_matrix = [quality_R1_P1, quality_R1_P2, quality_R2_P1, quality_R2_P2]

            print("%s" % Quality_matrix[select[0]], file=file_handle[0])
            print("%s" % Quality_matrix[select[1]], file=file_handle[1])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-m", help="minimum length", type=int, default=20)
    parser.add_argument("-n", help="match length", type=int, default=10)
    parser.add_argument("-prefix", default="Test", help="prefix")
    parser.add_argument("-m1", help="merge P1")
    parser.add_argument("-m2", help="merge P2")
    parser.add_argument("-A1", help="unmerge A1")
    parser.add_argument("-G1", help="unmerge G1")
    parser.add_argument("-A2", help="unmerge A2")
    parser.add_argument("-G2", help="unmerge G2")

    args = parser.parse_args()

    match_length = args.n
    min_length = args.m
    f_prefix = args.prefix

    linker_seq = "N" * match_length

    Out_ok1 = open(f_prefix + "_OK_R1.fastq", "w")
    Out_ok2 = open(f_prefix + "_OK_R2.fastq", "w")
    Out_short1 = open(f_prefix + "_Short_R1.fastq", "w")
    Out_short2 = open(f_prefix + "_Short_R2.fastq", "w")
    Out_discard1 = open(f_prefix + "_Discard_R1.fastq", "w")
    Out_discard2 = open(f_prefix + "_Discard_R2.fastq", "w")

    T1 = gzip.open(args.m1, "rb")
    T2 = gzip.open(args.m2, "rb")
    A1 = gzip.open(args.A1, "rb")
    G1 = gzip.open(args.G1, "rb")
    A2 = gzip.open(args.A2, "rb")
    G2 = gzip.open(args.G2, "rb")

    t1 = [i.rstrip() for i in T1.readlines()]
    t2 = [i.rstrip() for i in T2.readlines()]
    y1 = [i.rstrip() for i in G1.readlines()]
    y2 = [i.rstrip() for i in A1.readlines()]
    y3 = [i.rstrip() for i in G2.readlines()]
    y4 = [i.rstrip() for i in A2.readlines()]

    parser_no_merge_reads_deal(Out_ok1, Out_ok2, Out_short1, Out_short2, Out_discard1, Out_discard2, y1, y2, y3, y4)
    parser_merged_reads_deal(Out_ok1, Out_ok2, Out_discard1, Out_discard2, t1, t2)

    Out_ok1.close()
    Out_ok2.close()
    Out_short1.close()
    Out_short2.close()
    Out_discard1.close()
    Out_discard2.close()
    T1.close()
    T2.close()
    A1.close()
    G1.close()
    A2.close()
    G2.close()
