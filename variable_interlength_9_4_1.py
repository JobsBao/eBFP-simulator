import math
import random
import cmath
import copy
import time
import numpy as np
from matplotlib import pyplot as plt
import decimal


class decimal_HPIEEE754_transform():
    def __init__(self, bytes, exponent_bits, fraction_bits):
        self.bytes = bytes
        self.e_bits = exponent_bits #指数bit
        self.f_bits = fraction_bits #

    def generate_zero(self):  # 生成值为1的浮动块（[0,0,0,0,0,0,0,1]）
        res = [0] * self.f_bits

        return res
    def generate_zero_exp(self):  # 生成值为1的浮动块（[0,0,0,0,0,0,0,1]）
        res = [0] * self.e_bits
        return res
    def d_to_b(self, decimal):#高精度十进制数转BFP；输入参数：高精度十进制数；输出参数：BFP
        signal_flag = 0
        if decimal < 0:  # 确定标志位
            decimal = -decimal
            signal_flag = 1
        exponent_block = []  # 计算浮动块
        block_float_point = []
        exponent_block.append(signal_flag)  # 浮动块首位为符号位！！
        if decimal == 0:
            # block_exponent=0
            # exponent_block=[0,0,0,0,0,0,0,0]
            block_float_point.append(self.generate_zero_exp())
            for i in range(self.bytes - 1):
                block_float_point.append(self.generate_zero())
            return block_float_point
        else:
            block_exponent = math.ceil((math.floor(math.log(decimal, 2)) + 1) / self.f_bits)  # 计算浮动位
        max_digit = block_exponent * self.f_bits  # 计算最大浮动位
        min_digit = max_digit - (self.bytes - 1) * self.f_bits  # 计算最小浮动位
        decimal += (2 ** (min_digit - 1))  # 四舍五入
        block_exponent = math.ceil((math.floor(math.log(decimal, 2)) + 1) / self.f_bits)  # 四舍五入后重新计算浮动位
        max_digit = block_exponent * self.f_bits  # 重新计算最大浮动位
        min_digit = max_digit - (self.bytes - 1) * self.f_bits  # 重新计算最小浮动位
        block_exponent_754 = block_exponent + (2 ** (self.e_bits - 2) - 1)  # 根据754规则对浮动位进行上浮,得到指数位
        if block_exponent_754 <= 0:
            #print('block_exponent_754: ', block_exponent_754)
            raise EOFError('underflow!!')
        underflow_flag = 0  # 用于判断underflow
        for i in range(self.e_bits - 2, -1, -1):  # 存储指数位
            exponent_block.append(int(block_exponent_754 // (2 ** i)))
            if exponent_block[-1] > 1:
                #print('block_exponent_754: ', block_exponent_754)
                raise EOFError('overflow!!')
            block_exponent_754 -= exponent_block[-1] * (2 ** i)
            if exponent_block[-1] == 1:
                underflow_flag = 1
        if underflow_flag == 0:
            raise EOFError('underflow!!')

        # 最终输出（指数块+浮动块）
        block_float_point.append(exponent_block)  # 先将指数块加入

        fix_point_list = []  # 浮动位存储器（暂时保存所有浮动位）
        for i in range(max_digit - 1, min_digit - 1, -1):  # 填充浮动位存储器
            fix_point_list.append(int(decimal // (2 ** i)))
            decimal -= fix_point_list[-1] * (2 ** i)
        for i in range(self.bytes - 1):  # 将浮动位存储器内的浮动位切片装入浮动块
            block_float_point.append(fix_point_list[i * self.f_bits: (i + 1) * self.f_bits])
        return block_float_point

    def complex_d_to_b(self, complex_decimal):
        complex_BFP = []
        real_part = complex_decimal.real
        image_part = complex_decimal.imag
        real_BFP = self.d_to_b(real_part)
        image_BFP = self.d_to_b(image_part)
        complex_BFP.append(real_BFP)
        complex_BFP.append(image_BFP)
        return complex_BFP

    def b_to_d(self, block_float_point):#BFP转低精度十进制数；输入参数：BFP；输出参数：低精度十进制数
        exponent_block = block_float_point[0]#提取指数块
        exponent_bits = exponent_block[1:]  #指数块去掉符号位
        block_exponent_754 = 0#计算指数位
        A = range(len(exponent_bits))
        for i in range(len(exponent_bits)):
            block_exponent_754 += exponent_bits[-(i + 1)] * (2 ** i)
        block_exponent = block_exponent_754 - (2 ** (self.e_bits - 2) - 1) #十进制的指数位 1
        decimal = 0
        for i in range(self.bytes - 1):
            current_decimal = 0
            for j in range(len(block_float_point[i + 1])):#计算浮动位
                current_decimal += block_float_point[i + 1][-(j + 1)] * (2 ** j)  #从左到由计算
            decimal += current_decimal * 2 ** ((block_exponent - i - 1) * self.f_bits)#浮动位乘指数位
        if exponent_block[0] == 1:#确定正负号
            decimal = -decimal
        return decimal

    def complex_b_to_d(self,complex_BFP):
        # print('re',complex_BFP[0])
        real_part = self.b_to_d(complex_BFP[0])
        image_part = self.b_to_d(complex_BFP[1])
        complex_result = real_part+1j*image_part
        return complex_result

    def cd_to_b(self, A):
        B = []
        for i in range(len(A)):
            B.append([])
            for j in range(len(A[0])):
                #B[i].append([])
                B[i].append(self.complex_d_to_b(A[i][j]))
        return B

    def rd_to_b(self, A):
        B = []
        for i in range(len(A)):
            B.append([])
            for j in range(len(A[0])):
                #B[i].append([])
                B[i].append(self.d_to_b(A[i][j]))
        return B

    def cb_to_d(self, A):
        B = []
        for i in range(len(A)):
            B.append([])
            for j in range(len(A[0])):
                # print('a',A[i][j])
                tmp=self.complex_b_to_d(A[i][j])
                B[-1].append(tmp)
        return B

    def rb_to_d(self, A):
        B = []
        for i in range(len(A)):
            B.append([])
            for j in range(len(A[0])):

                tmp=self.b_to_d(A[i][j])
                B[-1].append(tmp)
        return B

    def d_to_list(self, lst):
        res = []
        for i in range(len(lst)):
            res.append(self.d_to_b((lst[i])))
        return res

    def list_to_d(self, D):
        res = []
        for i in range(len(D)):
            res.append(self.b_to_d((D[i])))
        return res

    def complex_to_real(self, complex_M):
        res = []
        for i in range(len(complex_M)):
            res.append([])
            for j in range(len(complex_M[0])):
                res[i].append(complex_M[i][j][0])
        return res

    def complex_to_imag(self, complex_M):
        res = []
        for i in range(len(complex_M)):
            res.append([])
            for j in range(len(complex_M[0])):
                res[i].append(complex_M[i][j][1])
        return res

    def real_to_complex(self,real_M):
        res=[]
        for i in range(len(real_M)):
            res.append([])
            for j in range(len(real_M[0])):
                res[i].append([real_M[i][j],self.d_to_b(0.0)])
        return res
    def simple_inverse_signal(self, B):
        res = copy.deepcopy(B)
        res[0][0] = res[0][0] ^ 1
        return res

    def inverse_signal(self, B):
        res = copy.deepcopy(B)
        res[0][0][0] = res[0][0][0] ^ 1
        res[1][0][0] = res[1][0][0] ^ 1
        return res

    def simple_Matrix_inverse_signal(self,M):
        res=[]
        for i in range(len(M)):
            res.append([])
            for j in range(len(M[0])):
                res[i].append(self.simple_inverse_signal(M[i][j]))
        return res

    def complex_Matrix_Form(self, re, im):
        res = []
        for i in range(len(re)):
            res.append([])
            for j in range(len(re[0])):
                res[i].append([re[i][j], im[i][j]])
        return res


def generate_d_H_trans(high, e_bits, f_bits):#（256，8，8）
    d_H_trans = []
    d_H_trans.append([])
    for i in range(1, high + 1):
        tmp = decimal_HPIEEE754_transform(i, e_bits, f_bits)
        d_H_trans.append(tmp)
    return d_H_trans


class bnpc():
    def __init__(self, e_bits, f_bits):#（8，1）
        self.e_bits = e_bits
        self.f_bits = f_bits

    def generate_one(self):#生成值为1的浮动块（[0,0,0,0,0,0,0,1]）
        res = [0] * self.f_bits
        res[-1] = 1
        return res
    def generate_allone_exp(self):#生成值为1的浮动块（[0,0,0,0,0,0,0,1]）
        res = [1] * self.e_bits

        return res
    def generate_one_exp(self):#生成值为1的浮动块（[0,0,0,0,0,0,0,1]）
        res = [0] * self.e_bits
        res[-1] = 1
        return res

    def generate_exp_bias(self):
        res = [1] * self.e_bits
        res[0]=0
        res[1]=0
        return res

    def block_add(self, a, b): #定点block相加
        if len(a) != len(b):
            raise EOFError('different len!')
        l = len(a)          #block的长度
        buffer = [0] * l
        res = [0] * l
        c = 0
        for i in range(-1, -l - 1, -1):   #i= -1:1:-9
            if a[i] + b[i] + buffer[i] == 0:
                res[i] = 0
            elif a[i] + b[i] + buffer[i] == 1:
                res[i] = 1
            elif a[i] + b[i] + buffer[i] == 2:
                if i - 1 < -l:
                    c = 1
                else:
                    buffer[i - 1] = 1
                res[i] = 0
            elif a[i] + b[i] + buffer[i] == 3:
                if i - 1 < -l:
                    c = 1
                else:
                    buffer[i - 1] = 1
                res[i] = 1
            else:
                raise EOFError('bit add error!')
        return res, c

    def block_shift(self, block_fix_point, key, c_l):  # 对定点部分进行移位使得指数匹配
        zero_block = [0] * self.f_bits
        res = []
        for i in range(key):
            res.append(zero_block)
        res = res + block_fix_point
        if len(res) <= c_l - 1:
            for i in range(c_l - 1 - len(res)):
                res.append(zero_block)
            c = 0
        else:
            if res[c_l - 1][0] == 1:
                c = 1
            else:
                c = 0
            res = res[:c_l - 1]
        if len(res) != c_l - 1:
            raise EOFError('output length error!')
        return res, c

    def get_list_block_exponent(self, bfp_list):#加法
        block_exponent = []
        max_exponent = -100000000000
        for i in range(len(bfp_list)):
            block_exponent_754 = 0
            exponent_bits = bfp_list[i][0][1:]
            for i in range(len(exponent_bits)):
                block_exponent_754 += exponent_bits[-(i + 1)] * (2 ** i)
            block_exponent.append(block_exponent_754 - (2 ** (self.e_bits - 2) - 1))
            if block_exponent[-1] > max_exponent:
                max_exponent = block_exponent[-1]
        return block_exponent, max_exponent

    def get_block_exponent_754(self, decimal):  # 加法
        block_exponent_754 = decimal + (2 ** (self.e_bits - 2) - 1)
        if block_exponent_754 < 0:
            raise EOFError('underflow!!')
        exponent_byte = []
        for i in range(self.e_bits - 2, -1, -1):
            exponent_byte.append(int(block_exponent_754 // (2 ** i)))
            if exponent_byte[-1] > 1:
                raise EOFError('overflow!!')  # decimal过大, 无法保存
            block_exponent_754 -= exponent_byte[-1] * (2 ** i)
        return exponent_byte

    def block_compare(self, a, b):
        if len(a) != len(b):
            raise EOFError('different len!!!')
        for i in range(len(a)):
            if a[i] > b[i]:
                return 0
            if a[i] < b[i]:
                return 1
        return 2

    def compare(self, a, b):
        length = min(len(a), len(b))
        signal_a = a[0][0]
        signal_b = b[0][0]
        if signal_a < signal_b:
            return 0
        if signal_b < signal_a:
            return 1
        exp_com = self.block_compare(a[0][1:], b[0][1:])
        if exp_com == 0:
            if signal_a == 0:
                return 0
            else:
                return 1
        if exp_com == 1:
            if signal_a == 0:
                return 1
            else:
                return 0
        for i in range(1, length):
            frac_com = self.block_compare(a[i], b[i])
            if frac_com == 0:
                if signal_a == 0:
                    return 0
                else:
                    return 1
            if frac_com == 1:
                if signal_a == 0:
                    return 1
                else:
                    return 0
        if len(a) > len(b):
            if signal_a == 0:
                return 0
            else:
                return 1
        if len(a) < len(b):
            if signal_a == 0:
                return 1
            else:
                return 0
        return 2  ################################################################## 以上全部替换

    # 已替换
    def simple_add_same_signal(self, add_list, c_operate):#实数BFP相加；输入参数：add_list实数BFP列表(同号), 输出eBFP块数；输出参数：实数BFP
        # print("addlist",add_list)
        c_l = c_operate[1]
        res = []#输出BFP
        signal_flag = add_list[0][0][0]#确定符号位
        num = len(add_list)#确定待加数个数
        x_list = []
        for i in range(num):#提取每个待加数的块数
            x_list.append(len(add_list[i]))
        for i in range(num):#确保所有待加数符号位相同
            if add_list[i][0][0] != signal_flag:
                raise EOFError('error1!!')
        block_exponent, max_exponent = self.get_list_block_exponent(add_list)#提取各待加数指数位，得到最大指数位
        add_pool = [] #构建待加数存储池
        for i in range(c_l + max(num - 2, 0)): #add_pool[0]为进位池
            add_pool.append([])
        for i in range(num): #填充add_pool
            block_fix_point = add_list[i][1:]#提取所有浮动块
            shifted_block_fix_point, c = \
                self.block_shift(block_fix_point, max_exponent - block_exponent[i], c_l)#根据浮动位对浮动块移位
            for j in range(-1, -c_l, -1):#将浮动块对位放入add_pool
                add_pool[j].append(shifted_block_fix_point[j])
            if c == 1:#将进位放入add_pool
                add_pool[-1].append(self.generate_one())
        #print('add_pool', add_pool)
        for i in range(-1, -c_l - 1 - max(num - 2, 0), -1): #消化add_pool
            while len(add_pool[i]) > 1:#当对应位置消化至只有一个浮动块时停止消化
                if i > -(c_l - 1) and len(add_pool[i]) < 1:#浮动块个数为0
                    raise EOFError('error2!!')
                if i == -c_l and len(add_pool[i]) == 0:#进位池消化结束
                    break
                add_res, c = self.block_add(add_pool[i][-1], add_pool[i][-2])#浮动位相加
                del add_pool[i][-1]#删除已加数
                del add_pool[i][-1]
                add_pool[i].append(add_res)#添加加法结果
                if c == 1:
                    if i - 1 < -c_l - max(num - 2, 0):#不正常进位（进位池发生进位）
                        #print("add_list", add_list)
                        print(c_l)
                        raise EOFError('error3!!')
                    else:#添加进位
                        add_pool[i - 1].append(self.generate_one())
        c_bits = 0
        for i in range(-c_l, -c_l - 1 - max(num - 2, 0), - 1):
            if len(add_pool[i]) > 0:
                c_bits += 1
            else:
                break
        if c_bits > 0:
            if add_pool[(num - 1) - c_bits + c_l - 1][0][0] == 1:#四舍五入
                add_pool[(num - 1) - c_bits + c_l - 2].append(self.generate_one())
                for i in range(-1, -c_l - 1 - max(num - 2, 0), -1):  # 消化add_pool
                    while len(add_pool[i]) > 1:  # 当对应位置消化至只有一个浮动块时停止消化
                        if i > -(c_l - 1) and len(add_pool[i]) < 1:  # 浮动块个数为0
                            raise EOFError('error2!!')
                        if i == -c_l and len(add_pool[i]) == 0:  # 进位池消化结束
                            break
                        add_res, c = self.block_add(add_pool[i][-1], add_pool[i][-2])  # 浮动位相加
                        del add_pool[i][-1]  # 删除已加数
                        del add_pool[i][-1]
                        add_pool[i].append(add_res)  # 添加加法结果
                        if c == 1:
                            if i - 1 < -c_l - max(num - 2, 0):  # 不正常进位（进位池发生进位）
                                #print("add_list", add_list)
                                print(c_l)
                                raise EOFError('error3!!')
                            else:  # 添加进位
                                add_pool[i - 1].append(self.generate_one())
            c_bits = 0
            for i in range(-c_l, -c_l - 1 - max(num - 2, 0), - 1):
                if len(add_pool[i]) > 0:
                    c_bits += 1
                else:
                    break
        max_exponent += c_bits#存在进位，指数块加1
        exponent_byte = self.get_block_exponent_754(max_exponent)#转化为IEEE754格式
        exponent_byte.insert(0, signal_flag)#插入符号位
        block_fix_point = add_pool[(max(num - 2, 0) + 1) - c_bits: (max(num - 2, 0) + 1) - c_bits + c_l - 1]#提取所有浮动块
        res.append(exponent_byte)#输出结果添加指数块
        for i in range(len(block_fix_point)):#输出结果依次添加浮动块
            res.append(block_fix_point[i][0])
        if len(res) != c_l:
            raise EOFError('res length error!')
        return res

    def block_intershift(self, a, key):  # 8位bit,左移为key位，乘法
        alen = len(a)
        if key == 0:
            l = a
            h = [0] * self.f_bits
            return h, l
        elif key > alen:
            h = a[(key-alen):]+[0]*(key-alen)
            #t =
            l = [0] * self.f_bits
        else:
            h = [0]*(alen-key)+a[:key]
            l = a[key:]+[0]*key

            return h, l
    def block_time(self, a, b):  # 定点block相乘
        if len(a) != len(b):
            raise EOFError('different len!')
        l = len(a)  # block的长度
        low = [0] * l  # 低八位
        high = [0] * l  # 高八位    #初始化，不予考虑复杂度
        # 利用移位相加的方法
        for i in range(-1, -l - 1, -1):
            if b[i] == 1:  # 左移-i-1位
                t_high, t_low = self.block_intershift(a, -i - 1)  # eBFP: 可能的移位数：0,1,2,...,l-1 复杂度为：l(l-1)S/2,
                #  可能的赋值数：1,2,3,...,l-1, 复杂度：l(l-1)A/2
                low, c = self.block_add(low, t_low)   # block_add复杂度记为 eBFP: (l-1)*1.5L+L
                if c == 1:  # eBFP: 存在进位的可能性 0.5*2L = L
                    one = self.generate_one()
                    high, hc = self.block_add(high, one)
                high, hc = self.block_add(high, t_high)  # block_add复杂度记为 eBFP: (l-1)*1.5L+L
        # 循环平均次数：l/2
        # 每次循环复杂度：BFP：28S+28A+(11.5L+11.5L+L)
        # 每次循环复杂度 ：IEEE: l(l-1)A/2+l(l-1)S/2 +(l)*3L
        # 总复杂度:BFP：(28S+28A+24L)*l/2
        # IEEE 754: (l(l-1)A/2+l(l-1)S/2 +(l)*3L)*l/2
        mm = 1
        return high, low
    def block_minus(self, a, b): #定点block相减
        if len(a) != len(b):
            raise EOFError('different len!')
        l = len(a)
        buffer = [0] * l
        res = [0] * l # 最终结果
        brw = 0 #借位
        for i in range(-1, -l - 1, -1):
            if a[i] - b[i] - buffer[i] == 0:
                res[i] = 0
            elif a[i] - b[i] - buffer[i] == 1:
                res[i] = 1
            elif a[i] - b[i] - buffer[i] == -1:
                if i - 1 < -l:
                    brw = 1
                else:
                    buffer[i - 1] = 1
                res[i] = 1
            elif a[i] - b[i] - buffer[i] == -2:
                if i - 1 < -l:
                    brw = 1
                else:
                    buffer[i - 1] = 1
                res[i] = 0
            else:
                raise EOFError('bit add error!')
        return res, brw
    def simple_add(self, add_list, c_operate):#实数BFP相加；输入参数：实数BFP列表, 输出eBFP块数；输出参数：实数BFP
        '''for i in range(len(add_lists)):
            if self.b_to_d(add_lists[i]) != 0:
                add_list.append(add_lists[i])
        if len(add_list) == 0:
            return self.d_to_b(0)'''
        plus_signal = []
        minus_signal = []
        #print("------")
        for i in range(len(add_list)):
            if add_list[i][0][0] == 0:
                plus_signal.append(add_list[i])
            else:
                minus_signal.append(add_list[i])
        if len(plus_signal) == 0:
            return self.simple_add_same_signal(minus_signal, c_operate)
        elif len(minus_signal) == 0:
            return self.simple_add_same_signal(plus_signal, c_operate)
        else:
            plus_res = self.simple_add_same_signal(plus_signal, c_operate)
            minus_res = self.simple_add_same_signal(minus_signal, c_operate)
            if plus_res[0][0] == minus_res[0][0]:
                raise EOFError('same signal??')
            minus_res[0][0] = 0
            res = self.simple_minus([plus_res, minus_res],c_operate)
            for j in range(len(res) - 1):
                first_zero = 0
                for i in range(self.f_bits):
                    if res[1][i] != 0:
                        first_zero = 1
                        break
                if first_zero == 0:
                    del res[1]
                    res.append([0] * self.f_bits)
                    res[0], c = self.block_add(res[0], self.generate_one_exp())
                else:
                    break
            # for j in range(self.bytes - 1):
            #     first_zero = 0
            #     for i in range(self.bits):
            #         if res[1][i] != 0:
            #             first_zero = 1
            #             break
            #
            #     if first_zero == 0:
            #         del res[1]
            #         res.append([0, 0, 0, 0, 0, 0, 0, 0])
            #         res[0], c = self.block_add(res[0], [1, 1, 1, 1, 1, 1, 1, 1])
            #     else:
            #         break
            return res

    def simple_minus(self, minus_list, c_operate):  # 实数BFP相减；输入参数：实数BFP列表；输出参数：实数BFP
        c_l = c_operate[2]  # 减法操作位数
        e_len_a = len(minus_list[0][0])
        e_len_b = len(minus_list[1][0])
        len_a = len(minus_list[0])
        len_b = len(minus_list[1])
        f_len_a = len(minus_list[0][1])
        f_len_b = len(minus_list[1][1])
        # print("e_len_a",e_len_a)
        if decimal_HPIEEE754_transform(len_a, e_len_a, f_len_a).b_to_d(
                minus_list[0]) == 0.0 and decimal_HPIEEE754_transform(len_b, e_len_b, f_len_b).b_to_d(
                minus_list[1]) == 0.0:
            return decimal_HPIEEE754_transform(c_l, self.e_bits, self.f_bits).d_to_b(0)
        elif decimal_HPIEEE754_transform(len_b, e_len_b, f_len_b).b_to_d(minus_list[1]) == 0.0:
            dec = decimal_HPIEEE754_transform(len_a, e_len_a, f_len_a).b_to_d(minus_list[0])
            ans = decimal_HPIEEE754_transform(c_l, self.e_bits, self.f_bits).d_to_b(dec)

            return ans

        # print("minus_list[0]",minus_list[0])
        # print("minus_list[1]", minus_list[1])
        minus_list = copy.deepcopy(minus_list)
        res = []  # 输出BFP
        x_list = []
        signal_flag = minus_list[0][0][0]  # 确定符号位
        # print('signal_flag', signal_flag)
        num = len(minus_list)  # 确定待减数个数
        for i in range(num):  # 提取每个待加数的块数
            x_list.append(len(minus_list[i]))
        for i in range(num):  # 确保所有待加数符号位相同
            if minus_list[i][0][0] != signal_flag:
                minus_list[1][0][0] = minus_list[0][0][0]
                res = self.simple_add(minus_list,c_operate)
                return res
        minus_list_tmp = copy.deepcopy(minus_list)
        for i in range(len(minus_list_tmp)):
            minus_list_tmp[i][0][0] = 0
        if self.compare(minus_list_tmp[0], minus_list_tmp[1]) == 1:  # 比较被减数和减数的大小
            # print('d_H_trans.b_to_d(minus_list[0])', d_H_trans.b_to_d(minus_list[0]))
            # print('d_H_trans.b_to_d(minus_list[1])', d_H_trans.b_to_d(minus_list[1]))
            #print(d_H_trans[len(minus_list[0])].b_to_d(minus_list[0]))
            #print(d_H_trans[len(minus_list[1])].b_to_d(minus_list[1]))
            temp = minus_list[0]
            minus_list[0] = minus_list[1]
            minus_list[1] = temp
            # a = d_H_trans[x_list[1]].b_to_d(minus_list[0])
            # b = d_H_trans[x_list[0]].b_to_d(minus_list[1])
            # print('a', a)
            # print('b', b)
            # print('minus_list[0]', minus_list[0])
            # print('minus_list[1]', minus_list[1])
            if signal_flag == 0:
                new_flag = 1
            else:
                new_flag = 0
            signal_flag = new_flag  # 交换位置后符号位取反
            # print('signal_flag', signal_flag)
        block_exponent, max_exponent = self.get_list_block_exponent(minus_list)  # 提取各待减数指数位，得到最大指数位
        # print('max_exponent', max_exponent)
        minus_pool = []  # 构建待减数存储池
        for i in range(c_l):
            minus_pool.append([])
        # print('minus pool: ', minus_pool)
        c_list = []
        for i in range(num):  # 填充minus_pool
            block_fix_point = minus_list[i][1:]  # 提取所有浮动块
            shifted_block_fix_point, c = \
                self.block_shift(block_fix_point, max_exponent - block_exponent[i], c_l)  # 根据浮动位对浮动块移位
            c_list.append(c)
            if c == 1:  # 将四舍五入的进位加到对应的数字中
                '''temp1 = [[0 for i in range(self.bits)] for i in range(self.bytes)]  
                temp1[-1][-1] = 1  
                print('temp1', temp1)
                temp2 = temp1
                #print('shifted_block_fix_point', shifted_block_fix_point)
                for i in range(-1, -self.bytes, -1):
                    temp2[i] = shifted_block_fix_point[i]
                temp2[0] = [0 for i in range(self.bits)]
                print('temp2', temp2)
                rounding = bnp.simple_add([temp2, temp1])
                print('rounding', rounding)
                for i in range(-1, -self.bytes, -1):
                    shifted_block_fix_point[i] = rounding[i]
                print('shifted_block_fix_point', shifted_block_fix_point)'''
                temp1 = [0 for i in range(self.f_bits)]
                temp1[-1] = 1
                # print('temp1', temp1)
                shifted_block_fix_point[-1], c = self.block_add(shifted_block_fix_point[-1], temp1)
                # print('shifted_block_fix_point[-1]', shifted_block_fix_point[-1])
            for j in range(-1, -c_l, -1):  # 将浮动块对位放入minus_pool
                minus_pool[j].append(shifted_block_fix_point[j])

        minus_pool_copy = copy.deepcopy(minus_pool)

        for i in range(-1, -c_l, -1):
            '''print('i', i)
            print('minus_pool[i]', minus_pool[i])'''
            while len(minus_pool[i]) != 1:  # 当对应位置消化至只有一个浮动块时停止消化
                # print('len(minus_pool[i]): ', len(minus_pool[i]))
                if i > (-c_l) and len(minus_pool[i]) < 1:  # 浮动块个数为0
                    raise EOFError('error1!!')
                if i == (-c_l):  # minus_pool消化结束
                    break
                minus_res, brw = self.block_minus(minus_pool[i][0], minus_pool[i][1])  # 浮动位相减
                # print('minus_pool[i][0]', minus_pool[i][0])
                # print('minus_pool[i][1]', minus_pool[i][1])
                # print('minus_res', minus_res)
                # print('brw', brw)
                del minus_pool[i][1]  # 删除减数
                del minus_pool[i][0]
                minus_pool[i].append(minus_res)  # 添加减法结果
                if len(minus_pool[i]) > 1:
                    temp = minus_pool[i][0]
                    minus_pool[i][0] = minus_pool[i][1]
                    minus_pool[i][1] = temp
                # print('minus_pool[i]', minus_pool[i])
                if brw == 1:
                    if i < -c_l:  # 不正常借位
                        raise EOFError('error2!!')
                    else:  # 添加借位
                        minus_pool[i - 1].append(self.generate_one())
                        # print('minus_pool[i - 1]', minus_pool[i - 1])
        flag = 0
        if len(minus_pool[0]) > 0:  # 存在借位，更改指数块
            flag = 1
            # print('??')
            #print(minus_pool_copy)
            #print(c_list)
            #print(d_H_trans[len(minus_list[0])].b_to_d(minus_list[0]))
            #print(d_H_trans[len(minus_list[1])].b_to_d(minus_list[1]))
            # print('minus_pool[0]', minus_pool[0])
            max_exponent -= 1  # 存在借位，指数块减1
            '''print('max_exponent', max_exponent)
            print('signal_flag', signal_flag)'''
            exponent_byte = self.get_block_exponent_754(max_exponent)  # 转化为IEEE754格式
            exponent_byte.insert(0, signal_flag)  # 插入符号位
            block_fix_point = minus_pool[:-1]  # 提取所有浮动块
            res.append(exponent_byte)  # 输出结果添加指数块
            for i in range(len(block_fix_point)):  # 输出结果依次添加浮动块
                res.append(block_fix_point[i][0])
                # print('res', res)
        else:  # 不存在借位
            exponent_byte = self.get_block_exponent_754(max_exponent)
            exponent_byte.insert(0, signal_flag)
            block_fix_point = minus_pool[1:]
            res.append(exponent_byte)
            for i in range(len(block_fix_point)):
                res.append(block_fix_point[i][0])

        all_zero_flag = 0
        for j in range(len(res) - 1):
            first_zero = 0
            for i in range(self.f_bits):
                if res[1][i] != 0:
                    first_zero = 1
                    break

            if first_zero == 0:
                del res[1]
                res.append(self.generate_zero())
                res[0], c = self.block_add(res[0], self.generate_allone_exp())
                if j == len(res) - 2:
                    all_zero_flag = 1
            else:
                break
        if all_zero_flag == 1:
            res[0] = [0] * self.e_bits
        if flag == 1:
            # print('res', d_H_trans[len(res)].b_to_d(res))
            # print(res)
            pass
        return res

    def simple_time(self, a, b, c_operate):  # 实数BFP相加；输入参数：实数BFP列表, 输出eBFP块数；输出参数：实数BFP
        c_l = c_operate[3] # 最底层的乘法操作位数
        res = []  # 输出BFP
        signal_flag = a[0][0] ^ b[0][0]  # 确定乘法结果的符号位,IEEE754：L,eBFP：L
        x_list = []  # 确定乘数的块数的长度
        x_list.append(len(a))
        x_list.append(len(b))
        a_zero_input = []
        b_zero_input = []
        a_man = []
        b_man = []
        # a_zero_input.append([0] * len(a[0]))
        # b_zero_input.append([0] * len(b[0]))
        time_pool = []  # 构建乘法池
        # 对于输入存在0值的情况直接特殊处理
        for i in range(x_list[0] - 1):
            a_zero_input.append([0] * len(a[1]))
            a_man.append(a[i + 1])
        for i in range(x_list[1] - 1):
            b_zero_input.append([0] * len(b[1]))
            b_man.append(b[i + 1])
        if a_man == a_zero_input or b_man == b_zero_input:
            res.append([0] * self.e_bits)
            for i in range(0, c_l - 1, 1):
                t = [0] * self.f_bits
                res.append(t)
            zero = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)
            res = self.simple_add([res, zero], c_operate) ## 处理位数
            return res
        if c_l >= x_list[0] + x_list[1] - 1:
            for i in range(x_list[0] + x_list[1] - 2):  # add_pool[0]为进位池
                t = [0] * self.f_bits
                time_pool.append(t)
            # 计算尾数位 全精度保留计算结果
            # print('XXX:',x_list[0]+x_list[1])
            c = [0] * (x_list[0] + x_list[1])  # 进位
            one = self.generate_one()  # 初始化
            oc = [0] * (x_list[0] + x_list[1])  # 无法在当前循环中被处理的进位  #初始化部分不计入复杂度
            for j in range(-1, -x_list[1], -1):  # 乘数
                for i in range(-1, -x_list[0], -1):  # 被乘数
                    h, l = self.block_time(a[i], b[j])  # eBFP:112S+112A+96L
                    m1 = 0
                    m2 = 0
                    t1 = 0
                    t2 = 0
                    t3 = 0
                    z1 = 0
                    z2 = 0  # 标志位不计入复杂度
                    if oc[i + j + 1] == 1:  # eBFP: L
                        time_pool[i + j + 1], m2 = self.block_add(time_pool[i + j + 1], one)
                        oc[i + j + 1] = 0
                    time_pool[i + j + 1], m1 = self.block_add(time_pool[i + j + 1], l)  # m1，m2不会同时为1  #eBFP:11.5L
                    c[i + j] = m1 + m2  # 不考虑标志位
                    if oc[i + j] == 1:  # oc和c可能同时为1 # eBFP: L
                        time_pool[i + j], t3 = self.block_add(one, time_pool[i + j])
                        oc[i + j] = 0  # 处理完毕，进位值置为0
                    if c[i + j] == 1:  # eBFP: L
                        time_pool[i + j], t1 = self.block_add(one, time_pool[i + j])
                        c[i + j] = 0  # 处理完毕，进位值置为0, t1与t2 不能同时为0
                    time_pool[i + j], t2 = self.block_add(time_pool[i + j], h)  # t1与t2 不会同时为1 #eBFP:11.5L
                    c[i + j - 1] = t1 + t2  # eBFP: L
                    if c[i + j - 1] == 1:  # eBFP: L
                        time_pool[i + j - 1], z1 = self.block_add(time_pool[i + j - 1], one)
                        c[i + j - 1] = 0  # 处理完毕，进位值置为0
                    if t3 == 1:  # eBFP: L
                        time_pool[i + j - 1], z2 = self.block_add(time_pool[i + j - 1], one)
                        t3 = 0
                    oc[i + j - 2] = z1 + z2  # 不考虑标志位
            # 加上规格化判定，假如有尾数最高位有全零浮点块
            block_exponent, max_exponent = self.get_list_block_exponent([a, b])  # 提取各待乘数指数位，得到最大指数位（no need）

            # eBFP复杂度：7A+10L
            # IEEE 754复杂度：EA+(1.5E-0.5)L
            exponent = block_exponent[0] + block_exponent[1]  # eBFP:复杂度 7L #IEEE754:EL
            delnum = 0
            for i in range(c_l - 2, 0, -1):
                if time_pool[0] == [0] * self.f_bits:
                    del time_pool[0]
                    delnum = delnum + 1
                    exponent = exponent - 1
                else:
                    break
            # 不考虑复杂度
            # 加入四舍五入判定
            if c_l > (x_list[0] + x_list[1] - 2 - delnum):
                for i in range(0, c_l - x_list[1] - x_list[0] + 2 + delnum, 1):
                    t = [0] * self.f_bits
                    time_pool.append(t)

        elif c_l >= 1 & c_l < x_list[0] + x_list[1] - 1:
            for i in range(c_l + 1):  # add_pool[0]为进位池
                t = [0] * self.f_bits
                time_pool.append(t)
            # 计算尾数位 全精度保留计算结果
            # print('XXX:', c_l + 1)
            c = [0] * (c_l + 3)  # 进位
            one = self.generate_one()  # 初始化
            oc = [0] * (c_l + 3)  # 无法在当前循环中被处理的进位  #初始化部分不计入复杂度
            for j in range(-1, -x_list[1], -1):  # 乘数
                for i in range(-1, -x_list[0], -1):  # 被乘数
                    if i + j <= -(x_list[1] + x_list[0] - c_l - 1):
                        offset = x_list[1] + x_list[0] - c_l - 3
                        pos = i + j + offset
                        h, l = self.block_time(a[i], b[j])  # eBFP:112S+112A+96L
                        m1 = 0
                        m2 = 0
                        t1 = 0
                        t2 = 0
                        t3 = 0
                        z1 = 0
                        z2 = 0  # 标志位不计入复杂度
                        if oc[pos + 1] == 1:  # eBFP: L
                            time_pool[pos + 1], m2 = self.block_add(time_pool[pos + 1], one)
                            oc[pos + 1] = 0
                        time_pool[pos + 1], m1 = self.block_add(time_pool[pos + 1], l)  # m1，m2不会同时为1  #eBFP:11.5L
                        c[pos] = m1 + m2  # 不考虑标志位
                        if oc[pos] == 1:  # oc和c可能同时为1 # eBFP: L
                            time_pool[pos], t3 = self.block_add(one, time_pool[pos])
                            oc[pos] = 0  # 处理完毕，进位值置为0
                        if c[pos] == 1:  # eBFP: L
                            time_pool[pos], t1 = self.block_add(one, time_pool[pos])
                            c[pos] = 0  # 处理完毕，进位值置为0, t1与t2 不能同时为0
                        time_pool[pos], t2 = self.block_add(time_pool[pos], h)  # t1与t2 不会同时为1 #eBFP:11.5L
                        c[pos - 1] = t1 + t2  # eBFP: L
                        if c[pos - 1] == 1:  # eBFP: L
                            time_pool[pos - 1], z1 = self.block_add(time_pool[pos - 1], one)
                            c[pos - 1] = 0  # 处理完毕，进位值置为0
                        if t3 == 1:  # eBFP: L
                            time_pool[pos - 1], z2 = self.block_add(time_pool[pos - 1], one)
                            t3 = 0
                        oc[pos - 2] = z1 + z2  # 不考虑标志位

            # 加上规格化判定，假如有尾数最高位有全零浮点块
            block_exponent, max_exponent = self.get_list_block_exponent([a, b])  # 提取各待乘数指数位，得到最大指数位（no need）
            # print("block_exponent",block_exponent)
            # eBFP复杂度：7A+10L
            # IEEE 754复杂度：EA+(1.5E-0.5)L
            exponent = block_exponent[0] + block_exponent[1]  # eBFP:复杂度 7L #IEEE754:EL])
            delnum = 0
            for i in range(c_l - 2, 0, -1):
                if time_pool[0] == [0] * self.f_bits:
                    # print('time_pool1:', time_pool)
                    if (exponent < -2 ** (self.e_bits - 2) + 2):
                        break
                    del time_pool[0]
                    delnum = delnum + 1
                    exponent = exponent - 1

                else:
                    break
            # 不考虑复杂度
            # 加入四舍五入判定
            if c_l > (c_l + 1 - delnum):
                for i in range(0, c_l - c_l - 1 + delnum, 1):
                    t = [0] * self.f_bits
                    time_pool.append(t)

        else:

            print('error')
        temp_pool = time_pool
        if time_pool[c_l - 1][0] == 1 and time_pool[0][0] == 0:  # eBFP：0.5*2L = L #IEEE754:L
            for i in range(c_l - 2, -1, -1):
                time_pool[i], ms = self.block_add(time_pool[i], one)
                if ms == 0:
                    break
        exponent_byte = self.get_block_exponent_754(exponent)
        # print("exponent_byte", exponent_byte)
        # eBFP复杂度：7A+10L
        # IEEE 754复杂度：EA+(1.5E-0.5)L
        exponent_byte.insert(0, signal_flag)  # 添加符号位  #eBFP：1L #IEEE754:1L
        # 判断输出尾数是否为全0
        cl_man = []
        cl_zero_input = []
        for i in range(c_l - 1):
            cl_zero_input.append([0] * len(a[1]))
            cl_man.append(time_pool[i])
        if cl_zero_input == cl_man:
            res.append([0] * self.e_bits)
            for i in range(0, c_l - 1, 1):
                t = [0] * self.f_bits
                res.append(t)
            return res

        res.append(exponent_byte)  # 添加指数位  #eBFP复杂度：8A #IEEE754:EA

        # print('time_poor2',time_pool)
        for r in range(c_l - 1):
            res.append(time_pool[r])
        # print('res:', res)
        # 保留前byte-1个字节，添加到结果中  # eBFP复杂度：8(N-1)A #IEEE754:FA
        zero = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)
        res = self.simple_minus([res, zero], c_operate)
        return res
##################################################################### 以上全部替换
    def generate_zero(self):#除法
        res = [0] * self.f_bits
        return res

    def generate_zero_block(self,cl):#除法
        zero_block = []
        for i in range(cl-1):
            res=self.generate_zero()
            zero_block.append(res)
        return zero_block
    def in_block_left_shift(self,block,n):
        if n==0:
            return block
        else:
            res = [0] * self.f_bits
            res[0:-n]=block[n:]
            return res
    def in_block_right_shift(self,block,n):
        if n == 0:
            return block
        else:
            res = [0] * self.f_bits
            res[n:]=block[0:-n]
            return res
    def right_shift(self,bfp):
        for i in range(len(bfp)-1,-1,-1):
            if i!=0:
                bfp[i]=self.in_block_right_shift(bfp[i],1)
                add_b=self.in_block_left_shift(bfp[i-1],self.f_bits-1)
                # print(bfp[i],add_b)
                bfp[i],dd=self.block_add(bfp[i],add_b)
            else:
                bfp[i] = self.in_block_right_shift(bfp[i], 1)
        return bfp

    def left_shift(self,bfp):
        for i in range(len(bfp)):
            if i!=len(bfp)-1:
                bfp[i]=self.in_block_left_shift(bfp[i],1)
                add_b=self.in_block_right_shift(bfp[i+1],self.f_bits-1)
                bfp[i],dd=self.block_add(bfp[i],add_b)
            else:
                bfp[i] = self.in_block_left_shift(bfp[i], 1)
        return bfp

    def generate_one_block(self,l):#除法
        res_l = []
        for i in range(l):
            res = [0] * self.f_bits
            res_l.append(res)
        res_l[-1][-1] = 1
        return res_l

    def make_complement_code(self,D):  # 求补码，除法
        D_c = []
        for i in range(len(D)):
            if D[i] == 0:
                D_c.append(1)
            elif D[i] == 1:
                D_c.append(0)
        add_one = self.generate_one()
        complemnet_code, c = self.block_add(D_c, add_one)
        return complemnet_code

    def make_complement_code_exp(self,D):  # 求补码，除法
        D_c = []
        for i in range(len(D)):
            if D[i] == 0:
                D_c.append(1)
            elif D[i] == 1:
                D_c.append(0)
        add_one = self.generate_one_exp()
        # print('len',D)
        # print('D_C',D_c)
        # print('add_one', add_one)
        complemnet_code, c = self.block_add(D_c, add_one)
        return complemnet_code

    def make_block_complement_code(self,D):  # 求补码，除法
        D_c = []
        #print('D',D)
        for i in range(len(D)):
            D_c.append([])
            for j in range(self.f_bits):
                if D[i][j] == 0:
                    D_c[i].append(1)
                elif D[i][j] == 1:
                    D_c[i].append(0)
        add_one = self.generate_one_block(len(D))
        #print('add_one',add_one)
        complemnet_code = self.BFP_add(D_c, add_one)
        return complemnet_code

    def BFP_add(self,m,n):#除法   ####eBFP:84*L

        count=0
        f=0
        result=[]
        for i in range(-1,-len(m)-1,-1):
            res,c=self.block_add(m[i],n[i])
            if count==1:
                one_block=self.generate_one()
                res,f=self.block_add(one_block,res)
            count=max(c,f)
            result.insert(0,res)
        return result

    def qu_select_function(self,rSj):  # 基2——SRT算法商选择函数，除法
        #print("rSj",rSj)
        if self.f_bits==1:
            if rSj[0][0] == 0:
                if rSj[1][0] == 1 or rSj[2][0] == 1 or rSj[3][0] == 1:
                    qu = 1
                    # flag = '>=1/2'
                else:
                    qu = 0
                    # flag = '<1/2'
            elif rSj[0][0] == 1:
                if (rSj[1][0] == 0 or rSj[2][0] == 0):
                    qu = -1
                    # flag = '<-1/2'
                else:
                    # print('sds')
                    if rSj[3][0] == 1:
                        qu = 0
                        # flag = '>=-1/2'
                    else:
                        qu = -1
                        # flag = '<-1/2'
        elif self.f_bits==2:
            if rSj[0][0] == 0:
                if rSj[0][1] == 1 or rSj[1][0] == 1 or rSj[1][1] == 1:
                    qu = 1
                    #flag = '>=1/2'
                else:
                    qu = 0
                    #flag = '<1/2'
            elif rSj[0][0] == 1:
                if (rSj[0][1] == 0 or rSj[1][0] == 0):
                    qu = -1
                    #flag = '<-1/2'
                else:
                    #print('sds')
                    if rSj[1][1] == 1:
                        qu = 0
                        #flag = '>=-1/2'
                    else:
                        qu = -1
                        #flag = '<-1/2'
        elif self.f_bits==3:
            if rSj[0][0] == 0:
                if rSj[0][1] == 1 or rSj[0][2] == 1 or rSj[1][0] == 1:
                    qu = 1
                    #flag = '>=1/2'
                else:
                    qu = 0
                    #flag = '<1/2'
            elif rSj[0][0] == 1:
                if (rSj[0][1] == 0 or rSj[0][2] == 0):
                    qu = -1
                    #flag = '<-1/2'
                else:
                    #print('sds')
                    if rSj[1][0] == 1:
                        qu = 0
                        #flag = '>=-1/2'
                    else:
                        qu = -1
                        #flag = '<-1/2'
        else:
            if rSj[0][0] == 0:
                if rSj[0][1] == 1 or rSj[0][2] == 1 or rSj[0][3] == 1:
                    qu = 1
                    #flag = '>=1/2'
                else:
                    qu = 0
                    #flag = '<1/2'
            elif rSj[0][0] == 1:
                if (rSj[0][1] == 0 or rSj[0][2] == 0):
                    qu = -1
                    #flag = '<-1/2'
                else:
                    #print('sds')
                    if rSj[0][3] == 1:
                        qu = 0
                        #flag = '>=-1/2'
                    else:
                        qu = -1
                        #flag = '<-1/2'
        return qu
    def generate_one_position(self, pos):#除法
        res = [0] * self.f_bits
        res[pos] = 1
        return res
    def generate_one_block_position(self, byte_position, bit_position,cl):#除法
        res_l = []
        for i in range(cl-1):
            if i == byte_position:
                res = self.generate_one_position(bit_position)
            else:
                res = self.generate_zero()
            res_l.append(res)
        return res_l
    def in_block_arithmetic_right_shift_exp(self,block,n):
        flag=block[-1]
        res = [block[0]] * self.e_bits
        res[n:]=block[0:-n]
        return res,flag
    def BFP_add_sp(self, add_list, c_l):  # 实数BFP相加；输入参数：add_list实数BFP列表(同号), 输出eBFP块数；输出参数：实数BFP
        res = []  # 输出BFP
        signal_flag = add_list[0][0][0]  # 确定符号位
        num = len(add_list)  # 确定待加数个数
        x_list = []
        for i in range(num):  # 提取每个待加数的块数
            x_list.append(len(add_list[i]))
        for i in range(num):  # 确保所有待加数符号位相同
            if add_list[i][0][0] != signal_flag:
                raise EOFError('error1!!')
        block_exponent, max_exponent = self.get_list_block_exponent(add_list)  # 提取各待加数指数位，得到最大指数位
        add_pool = []  # 构建待加数存储池
        for i in range(c_l - 1):  # add_pool[0]为进位池
            add_pool.append([])
        for i in range(num):  # 填充add_pool
            block_fix_point = add_list[i][1:]  # 提取所有浮动块
            shifted_block_fix_point, c = \
                self.block_shift(block_fix_point, max_exponent - block_exponent[i], c_l)  # 根据浮动位对浮动块移位
            for j in range(-1, -c_l, -1):  # 将浮动块对位放入add_pool
                add_pool[j].append(shifted_block_fix_point[j])
            if c == 1:  # 将进位放入add_pool
                add_pool[-1].append(self.generate_one())
        i = -1  # 消化add_pool
        while True:
            while len(add_pool[i]) != 1:  # 当对应位置消化至只有一个浮动块时停止消化
                if i > -(c_l - 1) and len(add_pool[i]) < 1:  # 浮动块个数为0
                    raise EOFError('error2!!')
                if i == -c_l and len(add_pool[i]) == 0:  # 进位池消化结束
                    break
                add_res, c = self.block_add(add_pool[i][-1], add_pool[i][-2])  # 浮动位相加
                del add_pool[i][-1]  # 删除已加数
                del add_pool[i][-1]
                add_pool[i].append(add_res)  # 添加加法结果
                if c == 1:
                    if i - 1 < -len(add_pool):
                        pass
                    else:  # 添加进位
                        add_pool[i - 1].append(self.generate_one())
            i -= 1
            if i < -len(add_pool):
                break
        max_exponent += len(add_pool) - (c_l - 1)  # 存在进位，指数块加1
        exponent_byte = self.get_block_exponent_754(max_exponent)  # 转化为IEEE754格式
        exponent_byte.insert(0, signal_flag)  # 插入符号位
        block_fix_point = add_pool[:c_l - 1]  # 提取所有浮动块
        res.append(exponent_byte)  # 输出结果添加指数块
        for i in range(len(block_fix_point)):  # 输出结果依次添加浮动块
            res.append(block_fix_point[i][0])
        if len(res) != c_l:
            raise EOFError('res length error!')
        return res

    def simple_divide(self, a, b, c_operate):  # 实数除法——基2SRT算法
        cl = c_operate[4]  # 最底层除法操作位数
        q_list = self.generate_zero_block(cl)
        a_exponent_block = a[0]  # 提取指数块
        a_flag = a_exponent_block[0]
        a_exponent_bits = a_exponent_block[0:]  # 指数块去掉符号位（这里不需要去掉，因为最后都会被c_flag取代）
        b_exponent_block = b[0]  # 提取指数块
        b_flag = b_exponent_block[0]
        b_exponent_bits = b_exponent_block[0:]  # 指数块去掉符号位（这里不需要去掉，因为最后都会被c_flag取代）
        c_flag = a_flag ^ b_flag  ####eBFP:L*1 IEEE:L*1
        a_res = a[1:]
        b_res = b[1:]
        b_count = 0
        a_count = 0
        a_shift_flag = 0
        b_shift_flag = 0
        exponent_changed = 0  ####eBFP:8*N*A
        ####IEEE：（E+F)*A
        if self.f_bits == 1:
            if a_res[0][0] == 0:
                return decimal_HPIEEE754_transform(cl, self.e_bits, self.f_bits).d_to_b(0)
            if b_res[0][0] == 0:
                raise EOFError('被除数不能为0')

            a_res = self.right_shift(a_res)
            a_res = self.right_shift(a_res)
            a_res = self.right_shift(a_res)
            b_res = self.right_shift(b_res)
            b_res = self.right_shift(b_res)
            a_shift_flag = 3
            b_shift_flag = 2
        elif self.f_bits == 2:
            if a_res[0][0] == 0 and a_res[0][1] == 0:
                return decimal_HPIEEE754_transform(cl, self.e_bits, self.f_bits).d_to_b(0)
            if b_res[0][0] == 0 and b_res[0][1] == 0:
                raise EOFError('被除数不能为0')
            if a_res[0][0] == 0:
                a_res = self.right_shift(a_res)
                a_res = self.right_shift(a_res)
                a_shift_flag = 2
            elif a_res[0][0] != 0:
                a_res = self.right_shift(a_res)
                a_res = self.right_shift(a_res)
                a_res = self.right_shift(a_res)
                a_shift_flag = 3
            if b_res[0][0] == 0:
                b_res = self.right_shift(b_res)
                b_shift_flag = 1
            elif b_res[0][0] != 0:
                b_res = self.right_shift(b_res)
                b_res = self.right_shift(b_res)
                b_shift_flag = 2
        elif self.f_bits == 3:
            if a_res[0][0] == 0 and a_res[0][1] == 0 and a_res[0][2] == 0:
                return decimal_HPIEEE754_transform(cl, self.e_bits, self.f_bits).d_to_b(0)
            if b_res[0][0] == 0 and b_res[0][1] == 0 and b_res[0][2] == 0:
                raise EOFError('被除数不能为0')
            if a_res[0][0] == 0 and a_res[0][1] != 0:
                a_res = self.right_shift(a_res)
                a_res = self.right_shift(a_res)
                a_shift_flag = 2
            elif a_res[0][0] == 0 and a_res[0][1] == 0:
                a_res = self.right_shift(a_res)
                a_shift_flag = 1
            elif a_res[0][0] != 0:
                a_res = self.right_shift(a_res)
                a_res = self.right_shift(a_res)
                a_res = self.right_shift(a_res)
                a_shift_flag = 3
            if b_res[0][0] == 0 and b_res[0][1] != 0:
                b_res = self.right_shift(b_res)
                b_shift_flag = 1
            # elif b_res[0][0]==0 and b_res[0][1]==0:
            #     a_res = self.right_shift(a_res)
            #     a_shift_flag = 1
            elif b_res[0][0] != 0:
                b_res = self.right_shift(b_res)
                b_res = self.right_shift(b_res)
                b_shift_flag = 2
        else:
            i = 0
            while (a_res[0][i] != 1 and i <= self.f_bits - 2):
                a_count += 1
                i += 1
            if a_count == self.f_bits - 1 and a_res[0][self.f_bits - 1] == 0:
                return decimal_HPIEEE754_transform(cl, self.e_bits, self.f_bits).d_to_b(0)

            while (a_res[0][0] != 0 or a_res[0][1] != 0 or a_res[0][2] != 0):
                a_res = self.right_shift(
                    a_res)  ####eBFP:S*(8*(N-1)/7+2*8*(N-1)/7+3*8*(N-1)/7)=（48/7)*(N-1)*S   IEEE754:不需要此操作
                a_shift_flag += 1

            i = 0
            while (b_res[0][i] != 1 and i <= self.f_bits - 2):
                b_count += 1
                i += 1

            if b_count == self.f_bits - 1 and b_res[0][self.f_bits - 1] == 0:
                raise EOFError('被除数不能为0')

            if b_count < 2:
                b_shift_flag = 2 - b_count
                for k in range(b_shift_flag):
                    b_res = self.right_shift(b_res)
            elif b_count > 2:
                b_shift_flag = 2 - b_count
                for k in range(-b_shift_flag):
                    b_res = self.left_shift(
                        b_res)  ####eBFP:S(2*8*(N-1)/7+8*(N-1)/8+8*(N-1)/7+2*8*(N-1)/7+3*8*(N-1)/7+4*8*(N-1)/7+5*8*(N-1)/7)=(144/7)*(N-1)*S   IEEE754:不需要此操作
        # print('a_res',a_res)
        # print('b_res', b_res)
        b_complemnt = self.make_block_complement_code(b_res)  ####eBFP:L*(N-1)*8+2L
        ####IEEE:F*L+2L
        # print('b_com',b_complemnt)
        exponent_changed = a_shift_flag - b_shift_flag
        # print('ex',exponent_changed)
        for e in range(self.f_bits * (cl - 1)):
            a_res = self.left_shift(a_res)  #####eBFP:S*8*(N-1)
            #####IEEE:S*F
            q = self.qu_select_function(a_res)  # A

            # q_list.append(q)
            now_q = self.generate_one_block_position(e // self.f_bits, e % (self.f_bits), cl)
            if q == 1:
                a_res.insert(0, self.generate_one_exp())
                b_complemnt.insert(0, self.generate_one_exp())
                a_res = self.BFP_add_sp([a_res, b_complemnt], len(a_res))  ####eBFP:12*(N-1)*L
                ####IEEE:3/2*F*L
                a_res = a_res[1:]
                b_complemnt = b_complemnt[1:]
                q_list = self.BFP_add(q_list, now_q)  ####eBFP:12*(N-1)*L
                ####IEEE:3/2*F*L
            elif q == -1:
                ####eBFP:2/3*12*(N-1)*L+1/3*（L*(N-1)*8+2L)
                ####IEEE:2/3*3/2*F*L+1/3*(F*L+2L）
                a_res.insert(0, self.generate_one_exp())
                b_res.insert(0, self.generate_one_exp())
                # print('a_res before: ', a_res)
                # print('b_res before: ', b_res)
                a_res = self.BFP_add_sp([a_res, b_res], len(a_res))  ####eBFP:12*(N-1)*L
                a_res = a_res[1:]
                b_res = b_res[1:]
                # print('a_res after: ', a_res)
                # print('b_res after: ', b_res)
                ####IEEE:3/2*F*L
                now_q = self.make_block_complement_code(now_q)  ####eBFP:L*(N-1)*8+2L
                ####IEEE:F*L+2L
                q_list = self.BFP_add(q_list, now_q)  ####eBFP:12*(N-1)*L
                ####IEEE:3/2*F*L

        # print('q_l_1',q_list)
        jsq = 0
        if self.f_bits == 1:
            if (q_list[0][jsq] != 1):
                jsq += 1
        else:
            while (q_list[0][jsq] != 1):
                jsq += 1
            # print('canshu',jsq,exponent_changed)
        add_one = 0
        if exponent_changed > 0:
            if exponent_changed <= jsq:
                for s in range(exponent_changed):
                    q_list = self.left_shift(q_list)  ####eBFP:2*S*8*(N-1)
                    ####IEEE:不需要此操作
            else:
                add_one = 1
                first_block = self.generate_one()
                for r in range(exponent_changed):
                    first_block[-r - 1] = q_list[0][exponent_changed - 1 - r]
                # print('fs',first_block)
                for s in range(exponent_changed):
                    q_list = self.left_shift(q_list)
                jinwei = q_list[-1][0]
                del q_list[-1]
                q_list.insert(0, first_block)
                if jinwei == 1:
                    add_block = self.generate_one_block(cl - 1)
                    q_list = self.BFP_add(q_list, add_block)
        elif exponent_changed < 0:
            for s in range(-exponent_changed):
                q_list = self.right_shift(q_list)
        # print('q_l_2',q_list)
        b_exponent_bits_complement = self.make_complement_code_exp(b_exponent_bits)  ####eBFP:10*L
        ####IEEE:E*L+2*L
        c_exponent_bits, nos = self.block_add(a_exponent_bits, b_exponent_bits_complement)  ####eBFP:12*L
        ####IEEE:3/2*E*L
        c_exponent_bits, nos = self.block_add(c_exponent_bits, self.generate_exp_bias())  ####eBFP:12*L
        ####IEEE:3/2*E*L
        if add_one == 1:
            c_exponent_bits, nos = self.block_add(c_exponent_bits, self.generate_one_exp())

        c_exponent_bits[0] = c_flag

        q_list.insert(0, c_exponent_bits)  ####eBFP:8*A
        ####IEEE:E*A
        # print('c_ex', q_list)
        # return q_list,a_shift_flag,b_shift_flag,exponent_changed,c_exponent_bits
        for j in range(cl - 1):
            first_zero = 0
            for i in range(self.f_bits):
                if q_list[1][i] != 0:
                    first_zero = 1
                    break

            if first_zero == 0:
                del q_list[1]
                q_list.append(self.generate_zero())
                q_list[0], c = self.block_add(q_list[0], self.generate_one_exp())
            else:
                break
        return q_list

    def simple_sqrt_old(self, a, c_operate):  # 实数BFP开方；输入参数：实数BFP；输出参数：实数BFP
        c_l = c_operate[5]
        sq_list = self.generate_zero_block(cl)
        a_exponent_block = a[0]  # 提取指数块     ####eBFP:8*A IEEE:E*A
        # print('被开方数指数块',a_exponent_block)
        a_flag = a_exponent_block[0]  ###eBFP:1*A  IEEE:1*A
        if a_flag == 1:
            raise EOFError('被开方数不能是负数')
        a_exponent_bits = a_exponent_block[0:]  # 指数块去掉符号位（这里不需要去掉，因为最后都会被c_flag取代）
        bias_coplement = self.make_complement_code_exp(self.generate_exp_bias())  ###eBFP:10L IEEE:E*L+2L
        a_minus_bias, c = self.block_add(a_exponent_bits, bias_coplement)  ###eBFP:12L  EEE:1.5E*L
        # print('减去偏置指数块', a_minus_bias)
        a_minus_bias, flag = self.in_block_arithmetic_right_shift_exp(a_minus_bias, 1)  ####eBFP:8S IEEE:E*S

        a_minus_bias, c = self.block_add(a_minus_bias, self.generate_exp_bias())  ###eBFP:12L  EEE:1.5E*L
        if flag == 1:
            # print('进位')
            a_minus_bias, c = self.block_add(a_minus_bias, self.generate_one_exp())  ###eBFP:6L  EEE:0.75E*L
        a_res = a[1:]  ###eBFP:(N-1)*8L  IEEE:F*L
        # print('a_res',a_res)
        a_count = 0
        a_shift_num = 0
        # print('a_res_before_shift', a_res)
        while (a_res[0][a_count] != 1 and a_count <= self.f_bits - 2):
            a_count += 1
        if a_count == self.f_bits-1 and a_res[0][self.f_bits-1] == 0:
            return decimal_HPIEEE754_transform(cl, self.e_bits, self.f_bits).d_to_b(0)
        # print('a_count',a_count)
        # print('flag', flag)

        if self.f_bits==5 or 7:
            if flag==0:
                if a_count < 3:
                    while ((a_res[0][0] != 0) or (a_res[0][1] != 0) or (a_res[0][2] != 0) or (a_shift_num % 2 != 1)):
                        a_res = self.right_shift(a_res)  ###eBFP:0.5*(N-1)*8*S+0.5*(N-1)*8*3*S  IEEE:不需要此操作
                        a_shift_num = a_shift_num + 1
                elif a_count == 3:
                    a_res = self.right_shift(a_res)
                    a_shift_num = 1
                elif a_count == 4 or a_count == 5:
                    a_res = self.left_shift(a_res)
                    a_shift_num = -1
                elif a_count == 6 or a_count == 7:
                    a_res = self.left_shift(a_res)
                    a_res = self.left_shift(a_res)
                    a_res = self.left_shift(a_res)
                    a_shift_num = -3
            else:
                if a_count < 3:
                    while ((a_res[0][0] != 0) or (a_res[0][1] != 0) or (a_res[0][2] != 0) or (a_shift_num % 2 != 0)):
                        a_res = self.right_shift(a_res)  ###eBFP:0.5*(N-1)*8*S+0.5*(N-1)*8*3*S  IEEE:不需要此操作
                        a_shift_num = a_shift_num + 1
                elif a_count == 3 or a_count == 4:
                    # a_res = self.right_shift(a_res)
                    a_shift_num = 0
                elif a_count ==5 or 6:
                    a_res = self.left_shift(a_res)
                    a_res = self.left_shift(a_res)
                    a_shift_num = -2
                # elif a_count == 6 :
                #     a_res = self.left_shift(a_res)
                #     a_res = self.left_shift(a_res)
                #     a_res = self.left_shift(a_res)
                #     a_shift_num = -4

        else:
            if a_count < 3:
                while ((a_res[0][0] != 0) or (a_res[0][1] != 0) or (a_res[0][2] != 0) or (a_shift_num % 2 != 1)):
                    a_res = self.right_shift(a_res)  ###eBFP:0.5*(N-1)*8*S+0.5*(N-1)*8*3*S  IEEE:不需要此操作
                    a_shift_num = a_shift_num + 1
            elif a_count == 3:
                a_res = self.right_shift(a_res)
                a_shift_num = 1
            elif a_count == 4 or a_count == 5:
                a_res = self.left_shift(a_res)
                a_shift_num = -1
            elif a_count == 6 or a_count == 7:
                a_res = self.left_shift(a_res)
                a_res = self.left_shift(a_res)
                a_res = self.left_shift(a_res)
                a_shift_num = -3
        # print('a_res_after_shift', a_res)
        for i in range(self.f_bits * (cl - 1) - 3):
            # print('a_res', a_res)
            a_res = self.left_shift(a_res)  ###eBFP:(N-1)*8*S    IEEE:

            q = self.qu_select_function(a_res)  ###eBFP:8A       IEEE:8A
            if q == 1:
                old_sq_list = copy.deepcopy(sq_list)  ####     IEEE:
                sq_list[(i + 3) //self.f_bits][(i + 3) % self.f_bits] = 1  # q_list加1，直接改变对应位置为1即可      ####eBFP:     IEEE:
                # m = self.generate_one_block_position((i + 3) //self.f_bits, (i + 3) % self.f_bits,len(a_res)+1)  ####eBFP:     IEEE:
                m = self.generate_one_block_position((i + 3) // self.f_bits, (i + 3) % self.f_bits,
                                                     cl)  ####eBFP:     IEEE:
                if i != 0:
                    old_sq_list = self.left_shift(old_sq_list)  ###eBFP:    IEEE:
                    m = self.BFP_add(m, old_sq_list)  ###eBFP:    IEEE:
                dif=cl-1-len(a_res)
                if dif>0:
                    m=m[0:len(a_res)]
                elif dif<0:
                    for i in range(-dif):
                        m.append(self.generate_zero())
                m = self.make_block_complement_code(m)  ####eBFP:   ####IEEE:
                a_res = self.BFP_add(a_res, m)  ###eBFP:    IEEE:

            elif q == -1:  ####eBFP:2/3*((N-1)*8*A)+4/3*A+2/3*(N-1)*8*S+5/3*(N-1)*12L+L*(N-1)*8+2L
                # print('q',q)                                                                                                   ####ieee:2/3*F*A+4/3*A+2/3*F*S+5/3*1.5F*L+F*L+2L
                old_sq_list = copy.deepcopy(sq_list)  ####    IEEE:
                n_one = self.generate_one_block_position((i + 3) //self.f_bits, (i + 3) % self.f_bits,cl)  ####eBFP:     IEEE:
                n_one_com = self.make_block_complement_code(n_one)  ####eBFP:  IEEE:
                # print('n_one',n_one)
                sq_list = self.BFP_add(sq_list, n_one_com)  # 相当于q_list减去1，即与-1的补码相加  ###eBFP:    IEEE:
                m = self.generate_one_block_position((i + 3) //self.f_bits, (i + 3) % self.f_bits,cl)  ####eBFP:     IEEE:
                # print('m',m)
                m = self.make_block_complement_code(m)  ####eBFP:   ####IEEE:

                if i != 0:
                    old_sq_list = self.left_shift(old_sq_list)  ###eBFP:    IEEE:
                    m = self.BFP_add(m, old_sq_list)  ###eBFP:    IEEE:
                    # print('m',m)
                dif = cl - 1 - len(a_res)
                if dif > 0:
                    m = m[0:len(a_res)]
                elif dif < 0:
                    for i in range(-dif):
                        m.append(self.generate_zero())
                a_res = self.BFP_add(a_res, m)  ###eBFP:    IEEE:
            elif q == 0:
                sq_list[(i + 3) //  self.f_bits][(i + 3) %  self.f_bits] = 0

        # sq_list = self.left_shift(sq_list)
        # sq_list = self.left_shift(sq_list)
        # sq_list = self.left_shift(sq_list)
        #print('shift_num',a_shift_num)
        if self.f_bits==4:
            if a_shift_num == 1:#changed
                if flag == 0:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
            if a_shift_num == -1:
                if flag == 0:
                    sq_list = self.left_shift(sq_list)
                else:
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)

            if a_shift_num == 3:#changed
                if flag == 1:

                    sq_list = self.left_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
                else:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
            if a_shift_num == -3:
                # if flag==0:
                #     sq_list = self.right_shift(sq_list)
                #     sq_list = self.right_shift(sq_list)
                #     sq_list = self.right_shift(sq_list)
                if flag == 1:
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)


        elif self.f_bits==6:
            if a_shift_num == 1:#changed
                if flag == 0:
                    # sq_list = self.right_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                else:

                    sq_list = self.right_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
            if a_shift_num == -1:#chngde
                if flag == 0:
                    # sq_list = self.right_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                else:
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
            if a_shift_num == 3:#changde
                # if flag == 1:
                #     None
                #     # sq_list = self.right_shift(sq_list)
                #     # sq_list = self.right_shift(sq_list)
                #     # sq_list = self.right_shift(sq_list)
                #     # sq_list = self.right_shift(sq_list)
                if flag == 0:
                    # sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
            if a_shift_num == -3:
                # if flag==0:
                #     sq_list = self.right_shift(sq_list)
                #     sq_list = self.right_shift(sq_list)
                #     sq_list = self.right_shift(sq_list)
                if flag == 1:
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)


        if self.f_bits==5:
            if flag==0:
                if a_shift_num == 1:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                if a_shift_num==-1:
                    sq_list = self.left_shift(sq_list)
                if a_shift_num == 3:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
            if flag==1:
                if a_shift_num == 0:#changde
                    sq_list = self.right_shift(sq_list)
                # if a_shift_num == 2:#changde_right
                #     None
                if a_shift_num == 4:#changed_right
                    sq_list = self.left_shift(sq_list)


        if self.f_bits==7:
            if flag==0:
                if a_shift_num == 1:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                if a_shift_num==-1:
                    sq_list = self.left_shift(sq_list)
                if a_shift_num == 3:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
            if flag==1:
                if a_shift_num == -2:#changde_right
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                if a_shift_num == 0:#changde_right
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                if a_shift_num == 2:#changde_right
                    sq_list = self.right_shift(sq_list)
                if a_shift_num == 4:#changed_rigthy
                    # sq_list = self.left_shift(sq_list)
                    None



        '''
        if a_shift_num == 1:
            if flag == 0:
                # sq_list = self.right_shift(sq_list)
                sq_list = self.left_shift(sq_list)
                sq_list = self.left_shift(sq_list)
            else:
                # sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
                sq_list = self.right_shift(sq_list)
                sq_list = self.right_shift(sq_list)
        if a_shift_num == -1:
            if flag == 0:
                # sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
                sq_list = self.left_shift(sq_list)
            else:
                sq_list = self.right_shift(sq_list)
                sq_list = self.right_shift(sq_list)
                sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
        if a_shift_num == 3:
            if flag == 1:
                sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
            else:
                sq_list = self.left_shift(sq_list)
                sq_list = self.left_shift(sq_list)
                sq_list = self.left_shift(sq_list)
        if a_shift_num == -3:
            # if flag==0:
            #     sq_list = self.right_shift(sq_list)
            #     sq_list = self.right_shift(sq_list)
            #     sq_list = self.right_shift(sq_list)
            if flag == 1:
                sq_list = self.right_shift(sq_list)
                sq_list = self.right_shift(sq_list)
                sq_list = self.right_shift(sq_list)
                sq_list = self.right_shift(sq_list)
        '''
        sq_list.insert(0, a_minus_bias)
        return sq_list

    def simple_sqrt(self, a, c_operate):  # 实数BFP开方；输入参数：实数BFP；输出参数：实数BFP
        cl = c_operate[5]
        sq_list = self.generate_zero_block(cl)
        a_exponent_block = a[0]  # 提取指数块     ####eBFP:8*A IEEE:E*A
        # print('被开方数指数块',a_exponent_block)
        a_flag = a_exponent_block[0]  ###eBFP:1*A  IEEE:1*A
        if a_flag == 1:
            raise EOFError('被开方数不能是负数')
        a_exponent_bits = a_exponent_block[0:]  # 指数块去掉符号位（这里不需要去掉，因为最后都会被c_flag取代）
        bias_coplement = self.make_complement_code_exp(self.generate_exp_bias())  ###eBFP:10L IEEE:E*L+2L

        a_minus_bias, c = self.block_add(a_exponent_bits, bias_coplement)  ###eBFP:12L  EEE:1.5E*L
        # print('减去偏置指数块', a_minus_bias)
        a_minus_bias, flag = self.in_block_arithmetic_right_shift_exp(a_minus_bias, 1)  ####eBFP:8S IEEE:E*S

        a_minus_bias, c = self.block_add(a_minus_bias, self.generate_exp_bias())  ###eBFP:12L  EEE:1.5E*L
        if flag == 1:
            # print('进位')
            a_minus_bias, c = self.block_add(a_minus_bias, self.generate_one_exp())  ###eBFP:6L  EEE:0.75E*L
        a_res = a[1:]  ###eBFP:(N-1)*8L  IEEE:F*L
        # print('flag', flag)
        '''
        下面这一段是否必要？
        '''
        # dif = cl - 1 - len(a_res)
        # if dif > 0:
        #     for i in range(dif):
        #         a_res.append(self.generate_zero())
        # print('a_res',a_res)
        a_count = 0
        a_shift_num = 0
        if self.f_bits==1:
            if a_res[0][0]==0:
                return decimal_HPIEEE754_transform(cl, self.e_bits, self.f_bits).d_to_b(0)
            else:
                if flag==0:
                    a_res = self.right_shift(a_res)
                    a_res = self.right_shift(a_res)
                    a_res = self.right_shift(a_res)
                    a_shift_num = 3
                else:
                    a_res = self.right_shift(a_res)
                    a_res = self.right_shift(a_res)
                    a_res = self.right_shift(a_res)
                    a_res = self.right_shift(a_res)
                    a_shift_num = 4
        # print('a_res_before_shift', a_res)
        elif self.f_bits==2:
            if a_res[0][0]==0 and a_res[0][1]==0:
                return decimal_HPIEEE754_transform(cl, self.e_bits, self.f_bits).d_to_b(0)
            a_res = self.right_shift(a_res)
            a_res = self.right_shift(a_res)
            a_res = self.right_shift(a_res)
            a_shift_num = 3

        elif self.f_bits==3:
            if a_res[0][0]==0 and a_res[0][1]==0 and a_res[0][2]==0:
                return decimal_HPIEEE754_transform(cl, self.e_bits, self.f_bits).d_to_b(0)
            if flag==0:
                a_res = self.right_shift(a_res)
                a_res = self.right_shift(a_res)
                a_res = self.right_shift(a_res)
                a_shift_num = 3
            else:
                if a_res[0][0]==0 and a_res[0][1]==0:
                    a_res = self.right_shift(a_res)
                    a_res = self.right_shift(a_res)
                    a_shift_num = 2
                else:

                    a_res = self.right_shift(a_res)
                    a_res = self.right_shift(a_res)
                    a_res = self.right_shift(a_res)
                    a_res = self.right_shift(a_res)
                    a_shift_num = 4
        else:
            while (a_res[0][a_count] != 1 and a_count <= self.f_bits - 2):
                a_count += 1
            if a_count == self.f_bits - 1 and a_res[0][self.f_bits - 1] == 0:
                return decimal_HPIEEE754_transform(cl, self.e_bits, self.f_bits).d_to_b(0)
            #print('a_count', a_count)


            if self.f_bits == 5 or self.f_bits == 7:
                if flag == 0:
                    if a_count < 3:
                        while ((a_res[0][0] != 0) or (a_res[0][1] != 0) or (a_res[0][2] != 0) or (a_shift_num % 2 != 1)):
                            a_res = self.right_shift(a_res)  ###eBFP:0.5*(N-1)*8*S+0.5*(N-1)*8*3*S  IEEE:不需要此操作
                            a_shift_num = a_shift_num + 1
                    elif a_count == 3:
                        a_res = self.right_shift(a_res)
                        a_shift_num = 1
                    elif a_count == 4 or a_count == 5:
                        a_res = self.left_shift(a_res)
                        a_shift_num = -1
                    elif a_count == 6 or a_count == 7:
                        a_res = self.left_shift(a_res)
                        a_res = self.left_shift(a_res)
                        a_res = self.left_shift(a_res)
                        a_shift_num = -3
                else:
                    if a_count < 3:
                        while ((a_res[0][0] != 0) or (a_res[0][1] != 0) or (a_res[0][2] != 0) or (a_shift_num % 2 != 0)):
                            a_res = self.right_shift(a_res)  ###eBFP:0.5*(N-1)*8*S+0.5*(N-1)*8*3*S  IEEE:不需要此操作
                            a_shift_num = a_shift_num + 1
                    elif a_count == 3 or a_count == 4:
                        # a_res = self.right_shift(a_res)
                        a_shift_num = 0
                    elif a_count == 5 or 6:
                        a_res = self.left_shift(a_res)
                        a_res = self.left_shift(a_res)
                        a_shift_num = -2
                    # elif a_count == 6 :
                    #     a_res = self.left_shift(a_res)
                    #     a_res = self.left_shift(a_res)
                    #     a_res = self.left_shift(a_res)
                    #     a_shift_num = -4

            else:
                if a_count < 3:
                    while ((a_res[0][0] != 0) or (a_res[0][1] != 0) or (a_res[0][2] != 0) or (a_shift_num % 2 != 1)):
                        a_res = self.right_shift(a_res)  ###eBFP:0.5*(N-1)*8*S+0.5*(N-1)*8*3*S  IEEE:不需要此操作
                        a_shift_num = a_shift_num + 1
                elif a_count == 3:
                    a_res = self.right_shift(a_res)
                    a_shift_num = 1
                elif a_count == 4 or a_count == 5:
                    a_res = self.left_shift(a_res)
                    a_shift_num = -1
                elif a_count == 6 or a_count == 7:
                    a_res = self.left_shift(a_res)
                    a_res = self.left_shift(a_res)
                    a_res = self.left_shift(a_res)
                    a_shift_num = -3
        #print('a_res_after_shift', a_res)
        for i in range(self.f_bits * (cl - 1) - 3):
            # print('a_res', a_res)
            a_res = self.left_shift(a_res)  ###eBFP:(N-1)*8*S    IEEE:

            q = self.qu_select_function(a_res)  ###eBFP:8A       IEEE:8A

            # q_list.append(q)
            # print('a_res*2', a_res)
            if q == 1:
                old_sq_list = copy.deepcopy(sq_list)  ####     IEEE:
                sq_list[(i + 3) // self.f_bits][
                    (i + 3) % self.f_bits] = 1  # q_list加1，直接改变对应位置为1即可      ####eBFP:     IEEE:
                # m = self.generate_one_block_position((i + 3) //self.f_bits, (i + 3) % self.f_bits,len(a_res)+1)  ####eBFP:     IEEE:
                m = self.generate_one_block_position((i + 3) // self.f_bits, (i + 3) % self.f_bits,
                                                     cl)  ####eBFP:     IEEE:
                if i != 0:
                    old_sq_list = self.left_shift(old_sq_list)  ###eBFP:    IEEE:
                    m = self.BFP_add(m, old_sq_list)  ###eBFP:    IEEE:
                dif = cl - 1 - len(a_res)
                if dif > 0:
                    m = m[0:len(a_res)]
                elif dif < 0:
                    for i in range(-dif):
                        m.append(self.generate_zero())
                m = self.make_block_complement_code(m)  ####eBFP:   ####IEEE:
                a_res = self.BFP_add(a_res, m)  ###eBFP:    IEEE:

            elif q == -1:  ####eBFP:2/3*((N-1)*8*A)+4/3*A+2/3*(N-1)*8*S+5/3*(N-1)*12L+L*(N-1)*8+2L
                # print('q',q)                                                                                                   ####ieee:2/3*F*A+4/3*A+2/3*F*S+5/3*1.5F*L+F*L+2L
                old_sq_list = copy.deepcopy(sq_list)  ####    IEEE:
                n_one = self.generate_one_block_position((i + 3) // self.f_bits, (i + 3) % self.f_bits,
                                                         cl)  ####eBFP:     IEEE:
                n_one_com = self.make_block_complement_code(n_one)  ####eBFP:  IEEE:
                # print('n_one',n_one)
                sq_list = self.BFP_add(sq_list, n_one_com)  # 相当于q_list减去1，即与-1的补码相加  ###eBFP:    IEEE:
                m = self.generate_one_block_position((i + 3) // self.f_bits, (i + 3) % self.f_bits,
                                                     cl)  ####eBFP:     IEEE:
                # print('m',m)
                m = self.make_block_complement_code(m)  ####eBFP:   ####IEEE:

                if i != 0:
                    old_sq_list = self.left_shift(old_sq_list)  ###eBFP:    IEEE:
                    m = self.BFP_add(m, old_sq_list)  ###eBFP:    IEEE:
                    # print('m',m)
                dif = cl - 1 - len(a_res)
                if dif > 0:
                    m = m[0:len(a_res)]
                elif dif < 0:
                    for i in range(-dif):
                        m.append(self.generate_zero())
                a_res = self.BFP_add(a_res, m)  ###eBFP:    IEEE:
            elif q == 0:
                sq_list[(i + 3) // self.f_bits][(i + 3) % self.f_bits] = 0

        # sq_list = self.left_shift(sq_list)
        # sq_list = self.left_shift(sq_list)
        # sq_list = self.left_shift(sq_list)
        # print('shift_num', a_shift_num)
        if self.f_bits == 4:
            if a_shift_num == 1:  # changed
                if flag == 0:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
            if a_shift_num == -1:
                if flag == 0:
                    sq_list = self.left_shift(sq_list)
                else:
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)

            if a_shift_num == 3:  # changed
                if flag == 1:

                    sq_list = self.left_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
                else:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
            if a_shift_num == -3:
                # if flag==0:
                #     sq_list = self.right_shift(sq_list)
                #     sq_list = self.right_shift(sq_list)
                #     sq_list = self.right_shift(sq_list)
                if flag == 1:
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)


        elif self.f_bits == 6:
            if a_shift_num == 1:  # changed
                if flag == 0:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                else:

                    sq_list = self.right_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
            if a_shift_num == -1:  # chngde
                if flag == 0:

                    sq_list = self.left_shift(sq_list)
                else:
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)

            if a_shift_num == 3:  # changde
                # if flag == 1:

                if flag == 0:
                    # sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
            if a_shift_num == -3:
                # if flag==0:

                if flag == 1:
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)

        elif self.f_bits == 2:

            if a_shift_num == 3:  # 对的对的对的对的对的
                if flag == 1:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                if flag == 0:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)


        elif self.f_bits == 5:
            if flag == 0:
                if a_shift_num == 1:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                if a_shift_num == -1:
                    sq_list = self.left_shift(sq_list)
                if a_shift_num == 3:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
            if flag == 1:
                if a_shift_num == 0:  # changde
                    sq_list = self.right_shift(sq_list)
                # if a_shift_num == 2:#changde_right
                #     None
                if a_shift_num == 4:  # changed_right
                    sq_list = self.left_shift(sq_list)

        elif self.f_bits == 7:
            if flag == 0:
                if a_shift_num == 1:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                if a_shift_num == -1:
                    sq_list = self.left_shift(sq_list)
                if a_shift_num == 3:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
            if flag == 1:
                if a_shift_num == -2:  # changde_right
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                if a_shift_num == 0:  # changde_right
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                if a_shift_num == 2:  # changde_right
                    sq_list = self.right_shift(sq_list)
                if a_shift_num == 4:  # changed_rigthy
                    # sq_list = self.left_shift(sq_list)
                    None
        elif self.f_bits == 3:
            if a_shift_num==4:
                sq_list = self.left_shift(sq_list)
                sq_list = self.left_shift(sq_list)
            if a_shift_num==2:

                sq_list = self.left_shift(sq_list)
            if a_shift_num == 3:
                sq_list = self.left_shift(sq_list)
                sq_list = self.left_shift(sq_list)
                sq_list = self.left_shift(sq_list)
        elif self.f_bits == 1:
            if flag == 0:
                if a_shift_num == 1:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                if a_shift_num == -1:
                    sq_list = self.left_shift(sq_list)
                if a_shift_num == 3:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
            if flag == 1:
                if a_shift_num == 0:  # changde
                    sq_list = self.right_shift(sq_list)
                # if a_shift_num == 2:#changde_right
                #     None
                if a_shift_num == 4:  # changed_right
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
        elif self.f_bits==8:
            if a_shift_num == 1:
                if flag == 0:
                    # sq_list = self.right_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                else:
                    # sq_list = self.right_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
            if a_shift_num == -1:
                if flag == 0:
                    # sq_list = self.right_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                else:
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
            if a_shift_num == 3:
                if flag == 1:
                    sq_list = self.right_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
                    # sq_list = self.right_shift(sq_list)
                else:
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
                    sq_list = self.left_shift(sq_list)
            if a_shift_num == -3:
                # if flag==0:
                #     sq_list = self.right_shift(sq_list)
                #     sq_list = self.right_shift(sq_list)
                #     sq_list = self.right_shift(sq_list)
                if flag == 1:
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
                    sq_list = self.right_shift(sq_list)
        '''
        if a_shift_num == 1:
            if flag == 0:
                # sq_list = self.right_shift(sq_list)
                sq_list = self.left_shift(sq_list)
                sq_list = self.left_shift(sq_list)
            else:
                # sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
                sq_list = self.right_shift(sq_list)
                sq_list = self.right_shift(sq_list)
        if a_shift_num == -1:
            if flag == 0:
                # sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
                sq_list = self.left_shift(sq_list)
            else:
                sq_list = self.right_shift(sq_list)
                sq_list = self.right_shift(sq_list)
                sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
        if a_shift_num == 3:
            if flag == 1:
                sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
                # sq_list = self.right_shift(sq_list)
            else:
                sq_list = self.left_shift(sq_list)
                sq_list = self.left_shift(sq_list)
                sq_list = self.left_shift(sq_list)
        if a_shift_num == -3:
            # if flag==0:
            #     sq_list = self.right_shift(sq_list)
            #     sq_list = self.right_shift(sq_list)
            #     sq_list = self.right_shift(sq_list)
            if flag == 1:
                sq_list = self.right_shift(sq_list)
                sq_list = self.right_shift(sq_list)
                sq_list = self.right_shift(sq_list)
                sq_list = self.right_shift(sq_list)

        '''
        sq_list.insert(0, a_minus_bias)
        zero = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)
        q_list = self.simple_add([sq_list, zero], c_operate)
        return sq_list

    def add(self, addlist, c_operate): #复数BFP加法；输入参数：复数BFP列表；输出参数：复数BFP
        real_list = []
        image_list = []
        for i in range(len(addlist)):
            real_list.append(addlist[i][0])
            image_list.append(addlist[i][1])

        real_result=self.simple_add(real_list,c_operate)
        image_result = self.simple_add(image_list,c_operate)
        res = [] #res=sigma addlist
        res.append(real_result)
        res.append(image_result)
        return res

    def minus(self, a, b, c_operate):  # 复数BFP减法；输入参数：复数BFP a，复数BFP b；输出参数：复数BFP
        res = []  # res=a-b
        real_part=self.simple_minus([a[0],b[0]],c_operate)
        imag_part=self.simple_minus([a[1],b[1]],c_operate)
        res.append(real_part)
        res.append(imag_part)
        return res

    def time(self, a, b, c_operate):  # 复数BFP乘法；输入参数：复数BFP a，复数BFP b；输出参数：复数BFP
        res = []  # res=a*b
        a_rb_r=self.simple_time(a[0], b[0], c_operate)
        a_ib_i=self.simple_time(a[1], b[1], c_operate)
        print('a_rb_r',decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(a_rb_r))
        print('a_ib_i',decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(a_ib_i))
        real_part=self.simple_minus([a_rb_r, a_ib_i], c_operate)
        print('real_part',decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(real_part))
        print("##########参与运算的BFP")
        print('a_rb_r',a_rb_r)
        print('a_ib_i',a_ib_i)
        print('real_part',real_part)
        a_rb_i = self.simple_time(a[0], b[1], c_operate)
        a_ib_r = self.simple_time(a[1], b[0], c_operate)
        print('a_rb_i',decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(a_rb_i))
        print('a_ib_r',decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(a_ib_r))
        imag_part=self.simple_add([a_rb_i,a_ib_r],c_operate)
        print('imag_part',decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(imag_part))
        print("##########参与运算的BFP")
        print('a_rb_i',a_rb_i)
        print('a_ib_r',a_ib_r)
        print('imag_part',imag_part)

        res.append(real_part)
        res.append(imag_part)

        return res
    def divide(self, a, b, c_operate):  # 复数BFP除法；输入参数：复数BFP a，复数BFP b；输出参数：复数BFP
        x2 = self.simple_time(b[0],b[0],c_operate)
        y2 = self.simple_time(b[1],b[1],c_operate)
        # print('y2',y2)
        x2y2 = self.simple_add([x2,y2],c_operate)
        first_part = self.simple_divide(b[0], x2y2, c_operate)
        _y=self.simple_time(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(-1), b[1], c_operate)
        second_part = self.simple_divide(_y, x2y2, c_operate)
        a1 = [first_part,second_part]
        res = self.time(a, a1, c_operate)
        return res


    def sqrt(self, a, c_operate):#复数BFP开方；输入参数：复数BFP；输出参数：复数BFP
        if decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(a[1])==0:
            zero = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)
            x = self.simple_minus([self.simple_sqrt(a[0],c_operate),zero],c_operate)
            return [x,decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)]
        else:
            x2 = self.simple_time(a[0], a[0], c_operate)
            y2 = self.simple_time(a[1], a[1], c_operate)
            x2y2 = self.simple_add([x2, y2], c_operate)
            sqrt_x2y2 = self.simple_sqrt(x2y2, c_operate)
            sqrt_x2y2_plus_x = self.simple_add([sqrt_x2y2, a[0]],c_operate)
            sqrt_x2y2_plus_x_half = self.simple_divide(sqrt_x2y2_plus_x, decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(2), c_operate)
            sqrt_sqrt_x2y2_plus_x_half = self.simple_sqrt(sqrt_x2y2_plus_x_half,c_operate)
            sqrt_x2y2_minus_x = self.simple_minus([sqrt_x2y2, a[0]],c_operate)
            sqrt_x2y2_minus_x_half = self.simple_divide(sqrt_x2y2_minus_x, decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(2),c_operate)  ##这里可能要优化
            sqrt_sqrt_x2y2_minus_x_half = self.simple_sqrt(sqrt_x2y2_minus_x_half,c_operate)
            print('sqrt_sqrt_x2y2_minus_x_half',sqrt_sqrt_x2y2_minus_x_half)
            print('sqrt_sqrt_x2y2_plus_x_half', sqrt_sqrt_x2y2_plus_x_half)
            return [sqrt_sqrt_x2y2_plus_x_half, sqrt_sqrt_x2y2_minus_x_half]

    def simple_abs(self, a):
        res = copy.deepcopy(a)
        res[0][0] = 0
        # print('a: ', I_to_d(self.time(a, self.conj(a, param), param)[0], param))
        return res

    def conj(self, a):
        b = []
        b.append(a[0])
        img = copy.deepcopy(a[1])
        img[0][0] = 1 ^ img[0][0]
        b.append(img)
        return b

    def abs(self, a, c_operate):
        # print('a: ', I_to_d(self.time(a, self.conj(a, param), param)[0], param))
        return self.simple_sqrt(self.time(a, self.conj(a),c_operate)[0], c_operate)

    def simple_power_two(self, a, c_operate):
        #print('sssssssss',a)
        return self.simple_time(a, a,c_operate)

    def a_time_A(self, a, A, c_operate):
        res=copy.deepcopy(A)
        for i in range(len(res)):
            for j in range(len(res[0])):
                res[i][j] = self.time(a, res[i][j], c_operate)
        return res
    def simple_A_add_B(self, A, B, c_operate):
        res=copy.deepcopy(A)
        for i in range(len(A)):
            for j in range(len(A[0])):
                res[i][j] = self.simple_add([A[i][j], B[i][j]],c_operate)
        return res
    def simple_A_minus_B(self, A, B,c_operate):
        res=copy.deepcopy(A)
        for i in range(len(A)):
            for j in range(len(A[0])):
                res[i][j] = self.simple_minus([A[i][j], B[i][j]],c_operate)
        return res
    def A_add_B(self, A, B, c_operate): #c_output 输出位，运算位，在矩阵加法中，这两者等价
        res = copy.deepcopy(A)
        for i in range(len(A)):
            for j in range(len(A[0])):
                res[i][j] = self.add([A[i][j], B[i][j]],c_operate)
        return res
    def A_minus_B(self, A, B, c_operate): #
        res = copy.deepcopy(A)
        for i in range(len(A)):
            for j in range(len(A[0])):
                res[i][j] = self.minus(A[i][j], B[i][j], c_operate)
        return res

    def trace(self, A, c_operate):
        if len(A) != len(A[0]):
            raise EOFError('different len')
        res = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_d_to_b(0 + 0j)
        for i in range(len(A)):
            res = self.add([res, A[i][i]], c_operate)
        return res

    def simple_dot(self, A, B, c_operate):
        C = []
        if len(A[0]) != len(B):
            raise EnvironmentError('len(A[0]) != len(B)')
        for i in range(len(A)):
            C.append([])
            for j in range(len(B[0])):
                tmp = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)
                # print(I_to_d(tmp, param))
                for k in range(len(B)):
                    tmp1 = self.simple_time(A[i][k], B[k][j],c_operate)
                    # print('tmp1',tmp1)
                    tmp = self.simple_add([tmp1, tmp],c_operate)
                    # print('tmp',tmp)
                C[i].append(tmp)
        return C

    def dot(self, A, B, c_operate):
        C = []

        if len(A[0]) != len(B):
            raise EnvironmentError('len(A[0]) != len(B)')
        for i in range(len(A)):
            C.append([])
            # print('c6', C)
            for j in range(len(B[0])):
                tmp = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_d_to_b(0+0j)
                #print(I_to_d(tmp, param))
                for k in range(len(B)):
                    tmp1 = self.time(A[i][k], B[k][j], c_operate)
                    tmp = self.add([tmp1, tmp],c_operate)
                    # print('tmp after add',tmp)
                    # print('----')
                C[i].append(tmp)
        return C

    def T(self, A):
        B = []
        for i in range(len(A[0])):
            B.append([])
            for j in range(len(A)):
                B[i].append(A[j][i])
        return B

    def conjt(self, A):
        B = []
        for i in range(len(A[0])):
            B.append([])
            for j in range(len(A)):
                #print('A[j][i]: ', A[j][i])
                B[i].append(self.conj(A[j][i]))
        return B

    def conj_m(self, A):
        B = []
        for i in range(len(A)):
            B.append([])
            for j in range(len(A[0])):
                #print('A[j][i]: ', A[j][i])
                B[i].append(self.conj(A[i][j]))
        return B

    def norm(self, A, c_operate):
        nm=decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0.0)
        for i in range(len(A)):
            for j in range(len(A[0])):
                nm=self.simple_add([nm,self.simple_power_two(A[i][j][0], c_operate),self.simple_power_two(A[i][j][1],c_operate)],c_operate)
        nm=self.simple_sqrt(nm, c_operate)
        return [nm,decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0.0)]

    def det(self, A, c_operate):
        # print(self.flag)
        if len(A) <= 0:
            return None
        if len(A) == 1:
            return A[0][0]
        else:
            s = [decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0), decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)]
            for i in range(len(A)):
                n = [[row[a] for a in range(len(A)) if a != i] for row in A[1:]]
                if i % 2 == 0:
                    s = self.add([s, self.time(A[0][i], self.det(n,c_operate),c_operate)],c_operate)
                else:
                    s = self.minus(s, self.time(A[0][i], self.det(n,c_operate),c_operate),c_operate)
        return s

    def sector(self, A, row, col):
        B = []
        for i in range(len(A)):
            if i < row[0] or i > row[1]:
                continue
            B.append([])
            for j in range(len(A[0])):
                if j < col[0] or j > col[1]:
                    continue
                B[-1].append(A[i][j])
        return B

    def simple_inv(self, A, c_operate):
        m, n = len(A), len(A[0])
        L = np.zeros((m, n))
        U = np.zeros((m, n))
        #print(self.e_bits,self.f_bits)
        L = decimal_HPIEEE754_transform(c_operate[1],self.e_bits,self.f_bits).rd_to_b(L)
        U = decimal_HPIEEE754_transform(c_operate[1],self.e_bits,self.f_bits).rd_to_b(U)
        # print('U', U)
        for i in range(n):
            U[0][i]=A[0][i]
        for i in range(m):
            L[i][0]=self.simple_divide(A[i][0],U[0][0],c_operate)

        for i in range(1,m):
            # 求U
            for j in range(i, n):
                s = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)
                for k in range(i):
                    s =self.simple_add([s , self.simple_time(L[i][k] , U[k][j],c_operate)],c_operate)
                U[i][j] =self.simple_minus([A[i][j], s],c_operate)
            # 求L
            L[i][i] = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(1)
            for j in range(i + 1, m):
                s = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)
                for k in range(i):
                    s = self.simple_add([s, self.simple_time(L[j][k] , U[k][i],c_operate)],c_operate)
                L[j][i] = self.simple_divide(self.simple_minus([A[j][i], s],c_operate) , U[i][i],c_operate)
        X = np.zeros((m, n))
        Y = np.eye(m)
        X=decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).rd_to_b(X)
        Y=decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).rd_to_b(Y)
        for i in range(n):
            for k in range(i + 1, m):
                tm=self.simple_dot(self.sector(L,[k,k],[0,n-1]), self.sector(Y,[0,m-1],[i,i]),c_operate)[0][0]

                Y[k][i] = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).simple_inverse_signal(tm)

        # 求X
        for i in range(n):
            for k in range(m - 1, -1, -1):
                X[k][i] = self.simple_divide(self.simple_minus([Y[k][i],
                                                                self.simple_dot(self.sector(U, [k, k], [0, n - 1]),
                                                                                self.sector(X, [0, m - 1], [i, i]), c_operate)[
                                                                    0][0]], c_operate), U[k][k], c_operate)
        return X

    def complex_inv(self,c,c_operate):
        a=decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_to_real(c)
        b=decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_to_imag(c)
        # print('a',np.array(self.rb_to_d(a)))
        comm = self.simple_A_add_B(a , self.simple_dot(b, self.simple_dot(self.simple_inv(a,c_operate), b,c_operate),c_operate),c_operate)
        # print('comm', np.array(self.rb_to_d(comm)))
        pre_realpart = comm
        real_part = self.simple_inv(pre_realpart,c_operate)
        # print('real_part', np.array(self.rb_to_d(real_part)))
        imag_part = self.simple_dot(self.simple_dot(self.simple_inv(a,c_operate),b,c_operate), self.simple_inv(comm,c_operate),c_operate)
        imag_part=decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).simple_Matrix_inverse_signal(imag_part)
        res = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_Matrix_Form(real_part, imag_part)
        return res

    def Neumann_series(self,A,NR,NT,num_term,c_operate):
        dim=len(A)
        Ik = np.identity(dim, dtype=complex)
        Ik_BFP=decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).cd_to_b(Ik)
        theta = 1 / (NR + NT) * np.identity(dim, dtype=complex)
        theta_BFP = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).cd_to_b(theta)
        temp2 = self.A_minus_B(Ik_BFP, self.dot(theta_BFP,A,c_operate),c_operate)
        z_inv_BFP=self.A_add_B(Ik_BFP,temp2,c_operate)
        for i in range(1, num_term):
            for j in range(2 ** i - 2 ** (i - 1)):
                temp2=self.dot(temp2,self.A_minus_B(Ik_BFP,self.dot(theta_BFP,A,c_operate),c_operate),c_operate)
            tmp = self.A_add_B(Ik_BFP, temp2, c_operate)
            z_inv_BFP=self.dot(z_inv_BFP, tmp, c_operate)
        res=self.dot(z_inv_BFP,theta_BFP,c_operate)
        return res

    def make_house_vec(self, test,c_operate): # 已替换
        x = copy.deepcopy(test)
        print("x",x)
        value_x = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).cb_to_d(x)
        print("x",value_x)
        #print("make_house_vec")
        #print("self.sector(x, [0, len(x) - 1], [0, 0])",self.sector(x, [0, len(x) - 1], [0, 0]))
        normx = self.dot(self.conjt(self.sector(x, [0, len(x) - 1], [0, 0])),
                         self.sector(x, [0, len(x) - 1], [0, 0]),c_operate)
        normx = normx[0][0]
        normx = self.sqrt(normx,c_operate)
        z = copy.deepcopy(x)
        z[0][0] = self.minus(z[0][0], normx,c_operate)
        ZHa = self.dot(self.conjt(z), x,c_operate)
        ZHa = ZHa[0][0]
        if decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_b_to_d(ZHa)!=0+0j:
            ZHa1 = self.divide(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_d_to_b(1+0j), ZHa,c_operate)
        else:
            ZHa1=decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_d_to_b(0+0j)
        return z, ZHa1,normx

    def house_bidiag_explicit_UV(self, A, c_operate): # 已替换
        m, n = len(A), len(A[0])
        assert m >= n
        U, Vt = np.identity(m, dtype='complex128'), np.identity(n, dtype='complex128')
        U = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).cd_to_b(U)
        Vt = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).cd_to_b(Vt)
        for col in range(n):

            v, beta ,nx= self.make_house_vec(self.sector(A, [col, len(A) - 1], [col, col]),c_operate)
            print("v",v)
            print("v",decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).cb_to_d(v))
            print("beta",decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_b_to_d(beta))
            if decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_b_to_d(beta) != 0 + 0j and abs(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_b_to_d(beta) )< 10 ** (10):
                test_1 = self.a_time_A(beta, v,c_operate)
                print('test1',decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).cb_to_d(test_1))
                test_3 = self.conjt(v)
                print('test3',decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).cb_to_d(test_3))
                test_2 = self.dot(test_1, self.conjt(v),c_operate)
                print('test2', decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).cb_to_d(test_2))
                tmp1 = self.A_minus_B(self.sector(A, [col, len(A) - 1], [col + 1, len(A[0]) - 1]),
                                      self.dot(self.dot(self.a_time_A(beta, v,c_operate), self.conjt(v),c_operate),
                                               self.sector(A, [col, len(A) - 1], [col + 1, len(A[0]) - 1]),c_operate),c_operate)
                i_count = 0
                for i in range(col, len(A)):
                    j_count = 0
                    for j in range(col+1, len(A[0])):
                        A[i][j] = tmp1[i_count][j_count]
                        j_count += 1
                    i_count += 1
                A[col][col] = nx
                for i in range(col+1, len(A)):
                    A[i][col] =decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_d_to_b(0)
                tmp2 = self.A_minus_B(self.sector(U, [0, len(U) - 1], [col, len(U[0]) - 1]),
                                      self.a_time_A(self.conj(beta), self.dot(
                                          self.dot(self.sector(U, [0, len(U) - 1], [col, len(U[0]) - 1]), v,c_operate),
                                          self.conjt(v),c_operate),c_operate),c_operate)
                i_count = 0
                for i in range(0, len(U)):
                    j_count = 0
                    for j in range(col , len(U[0])):
                        U[i][j] = tmp2[i_count][j_count]
                        j_count += 1
                    i_count += 1
            if col <= n - 2:
                v, beta,nx = self.make_house_vec(self.conjt(self.sector(A, [col, col], [col + 1, len(A[0])])),c_operate)
                if decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_b_to_d(beta) != 0 + 0j and abs(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_b_to_d(beta) )< 10 ** (10):
                    temp2 = self.A_minus_B(self.sector(A, [col+1, len(A) - 1], [col + 1, len(A[0]) - 1]), self.dot(
                        self.a_time_A(self.conj(beta),
                                      self.dot(self.sector(A, [col+1, len(A) - 1], [col + 1, len(A[0]) - 1]), v,c_operate),c_operate),
                        self.conjt(v),c_operate),c_operate)
                    k_count = 0
                    for k in range(col+1, len(A)):
                        l_count = 0
                        for l in range(col + 1, len(A[0])):
                            A[k][l] = temp2[k_count][l_count]
                            l_count += 1
                        k_count += 1
                    A[col][col+1] = nx
                    for i in range(col + 2, len(A[0])):
                        A[col][i] = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_d_to_b(0)
                    temp3 = self.A_minus_B(self.sector(Vt, [0, len(Vt) - 1], [col + 1, len(Vt[0]) - 1]), self.dot(
                        self.a_time_A(self.conj(beta),
                                      self.dot(self.sector(Vt, [0, len(Vt) - 1], [col + 1, len(Vt[0]) - 1]), v,c_operate),c_operate),
                        self.conjt(v),c_operate),c_operate)
                    k_count = 0
                    for k in range(0, len(Vt)):
                        l_count = 0
                        for l in range(col + 1, len(Vt[0])):
                            Vt[k][l] = temp3[k_count][l_count]
                            l_count += 1
                        k_count += 1
        return U, A, Vt

    def QR_D(self, A, c_operate): #检查完毕
        m, n = len(A), len(A[0])
        assert m >= n
        U = np.identity(m, dtype='complex128')
        U = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).cd_to_b(U)
        for col in range(n):
            v, beta, nx = self.make_house_vec(self.sector(A, [col, len(A) - 1], [col, col]), c_operate)
            if decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_b_to_d(beta) != 0 + 0j and abs(
                    decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_b_to_d(beta)) < 10 ** (10):
                tmp1 = self.A_minus_B(self.sector(A, [col, len(A) - 1], [col + 1, len(A[0]) - 1]),
                                      self.dot(self.dot(self.a_time_A(beta, v, c_operate), self.conjt(v), c_operate),
                                               self.sector(A, [col, len(A) - 1], [col + 1, len(A[0]) - 1]), c_operate), c_operate)
                i_count = 0
                for i in range(col, len(A)):
                    j_count = 0
                    for j in range(col + 1, len(A[0])):
                        A[i][j] = tmp1[i_count][j_count]
                        j_count += 1
                    i_count += 1
                A[col][col] = nx

                for i in range(col + 1, len(A)):
                    A[i][col] = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_d_to_b(0)
                tmp2 = self.A_minus_B(self.sector(U, [0, len(U) - 1], [col, len(U[0]) - 1]),
                                      self.a_time_A(self.conj(beta), self.dot(
                                          self.dot(self.sector(U, [0, len(U) - 1], [col, len(U[0]) - 1]), v, c_operate),
                                          self.conjt(v), c_operate), c_operate), c_operate)
                i_count = 0
                for i in range(0, len(U)):
                    j_count = 0
                    for j in range(col, len(U[0])):
                        U[i][j] = tmp2[i_count][j_count]
                        j_count += 1
                    i_count += 1
        return U, A

    def csr(self,x, y,c_operate): # 已替换
        if (decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(y) == 0):
            c = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(1)
            s = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)
            rr = x
        else:
            if (abs(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(y)) > abs(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(x))):

                tao = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).simple_inverse_signal(self.simple_divide(x , y, c_operate))
                s = self.simple_sqrt(self.simple_add([decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(1) , self.simple_power_two(tao,c_operate)],c_operate),c_operate)
                rr = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).simple_inverse_signal(self.simple_time(y , s,c_operate))
                s = self.simple_divide(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(1) , s,c_operate)
                c = self.simple_time(s , tao,c_operate)
            else:
                tao = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).simple_inverse_signal(self.simple_divide(y, x,c_operate))
                c = self.simple_sqrt(self.simple_add([decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(1), self.simple_power_two(tao,c_operate)],c_operate),c_operate)
                rr = self.simple_time(x, c,c_operate)
                c = self.simple_divide(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(1), c,c_operate)
                s = self.simple_time(c, tao,c_operate)
        return c, s, rr

    def updatecsvv(self, c, s, U, k,c_operate): # 已替换
        n = len(U)
        for i in range(n):
            t = U[i][k]
            U[i][k] = self.minus(self.time([c, decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)], t,c_operate), self.time([s, decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)], U[i][k + 1],c_operate),c_operate)
            U[i][k + 1] = self.add([self.time([s, decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)], t,c_operate), self.time([c, decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)], U[i][k + 1],c_operate)],c_operate)
        return U

    def updatecsvv_last(self, c1, s1, c2, s2, v1, v2,c_operate): # 已替换
        n = len(v1)
        for i in range(n):
            t = v1[i][0]
            v1[i][0] = self.add([self.time([c1, decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)], t, c_operate), self.time([s1, decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)], v2[i][0],c_operate)],c_operate)
            v2[i][0] = self.add([self.time([s2, decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)], t, c_operate), self.time([c2, decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)], v2[i][0],c_operate)],c_operate)
        return v1, v2
    def mul_22_submatrix(self,G, U1, U2,V,j,c_operate):  # 已替换
        g11 = G[0][0]
        g12 = G[0][1]
        g22 = G[1][1]
        a = self.simple_time(g11,g11,c_operate)
        b=self.simple_add([self.simple_time(g12,g12,c_operate),self.simple_time(g22,g22,c_operate)],c_operate)
        c=self.simple_time(g11,g12,c_operate)
        # compute the Jacobi rotation which diagonalizes[[a,c],[c,b]]
        theta=self.simple_divide(self.simple_minus([b,a],c_operate),self.simple_time(c,decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(2),c_operate),c_operate)
        t = self.simple_divide(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(
            np.sign(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(theta))), (self.simple_add(
            [self.simple_abs(theta),
             self.simple_sqrt(self.simple_add([decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(1), self.simple_time(theta, theta,c_operate)],c_operate),c_operate)],c_operate)),c_operate)
        cs = self.simple_divide(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(1), self.simple_sqrt(
            self.simple_add(
                [decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(1), self.simple_time(t, t,c_operate)],c_operate),c_operate),c_operate)
        sn = self.simple_time(cs , t,c_operate)
        # update columns i and i +1 of B
        cssn=np.zeros((2,2))
        cssn=decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).rd_to_b(cssn)
        cssn[0][0]=cs
        cssn[0][1] = sn
        cssn[1][0] =decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).simple_inverse_signal(sn)
        cssn[1][1] =cs
        G=self.simple_dot(G,cssn,c_operate)
        V=self.updatecsvv(cs,sn,V,j,c_operate)
        alpha = self.simple_sqrt(self.simple_add([self.simple_power_two(G[0][0],c_operate), self.simple_power_two(G[1][0],c_operate)],c_operate),c_operate)
        beta = self.simple_sqrt(self.simple_add([self.simple_power_two(G[0][1],c_operate), self.simple_power_two(G[1][1],c_operate)],c_operate),c_operate)
        c1 = self.simple_divide(G[0][0] , alpha,c_operate)
        c2 = self.simple_divide(G[1][1] , beta,c_operate)
        s1 = self.simple_divide(G[1][0] , alpha,c_operate)
        s2 = self.simple_divide(G[0][1] , beta,c_operate)
        # update rows i and i +1 of B
        c1s1s2c2=np.zeros((2,2))
        c1s1s2c2 = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).rd_to_b(c1s1s2c2)
        c1s1s2c2[0][0] = c1
        c1s1s2c2[0][1] = s1
        c1s1s2c2[1][0] = s2
        c1s1s2c2[1][1] = c2
        G = self.simple_dot(c1s1s2c2,G,c_operate)
        # update the matrix U of left singular vectors
        U1, U2 = self.updatecsvv_last(c1, s1, c2, s2, U1, U2,c_operate)
        # return G, U1, U2
        return G,U1,U2,V

    def estimation_lower_bound_minimum_singular_value1(self,B,c_operate):  # 已替换
        m=len(B)
        n=len(B[0])
        lamda = np.zeros(n)
        lamda=decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_list(lamda)
        lamda[n - 1] = self.simple_abs(B[n-1][n-1])
        for j in range(n - 2, -1, -1):
            lamda[j] =self.simple_time(self.simple_abs(B[j][j]),self.simple_divide(lamda[j+1],self.simple_add([lamda[j+1],self.simple_abs(B[j][j+1])],c_operate),c_operate),c_operate)
        u = np.zeros(n)
        u = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_list(u)
        u[0] = self.simple_abs(B[0][0])
        for j in range(n - 1):
            u[j + 1] = self.simple_time(self.simple_abs(B[j+1][j+1]),self.simple_divide(u[j],self.simple_add([u[j],self.simple_abs(B[j][j+1])],c_operate),c_operate),c_operate)
        B_Infinity = min(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).list_to_d(lamda))
        B_1 = min(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).list_to_d(u))
        sigma_low_bound = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(min(B_Infinity, B_1))
        return sigma_low_bound, lamda, u

    def QR_Wilkinson_shift_Iteration_once(self, tri, U, V, i_min, i_max, shift, c_operate): # 已替换
        n = i_max + 1
        x = self.simple_minus([self.simple_power_two(tri[i_min][i_min],c_operate), shift],c_operate)
        y = self.simple_time(tri[i_min][i_min], tri[i_min][i_min + 1],c_operate)
        for k in range(i_min, n - 1, 1):
            c, s, r = self.csr(x, y,c_operate)
            fir = np.zeros((2, 2))
            fir = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).rd_to_b(fir)
            fir[0][0] = tri[k][k]
            fir[0][1] = tri[k][k + 1]
            fir[1][0] = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)
            fir[1][1] = tri[k + 1][k + 1]
            sec = np.zeros((2, 2))
            sec = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).rd_to_b(sec)
            sec[0][0] = c
            sec[0][1] = s
            sec[1][0] = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).simple_inverse_signal(s)
            sec[1][1] = c
            rest = self.simple_dot(fir, sec,c_operate)
            x = rest[0][0]
            tri[k][k + 1] = rest[0][1]
            y = self.simple_time(fir[1][1], sec[1][0],c_operate)
            tri[k + 1][k + 1] = rest[1][1]
            V = self.updatecsvv(c, s, V, k,c_operate)
            if (k > 0):
                tri[k - 1][k] = r
            c, s, r = self.csr(x, y,c_operate)
            tri[k][k] = r
            U = self.updatecsvv(c, s, U, k,c_operate)
            if k != n - 2:
                fir = np.zeros((2, 2))
                fir = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).rd_to_b(fir)
                fir[0][0] = tri[k][k + 1]
                fir[0][1] = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)
                fir[1][0] = tri[k + 1][k + 1]
                fir[1][1] = tri[k + 1][k + 2]
                sec = np.zeros((2, 2))
                sec = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).rd_to_b(sec)
                sec[0][0] = c
                sec[0][1] = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).simple_inverse_signal(s)
                sec[1][0] = s
                sec[1][1] = c
                rest = self.simple_dot(sec, fir,c_operate)
                x = rest[0][0]
                y = self.simple_time(sec[0][1], fir[1][1],c_operate)
                tri[k + 1][k + 1] = rest[1][0]
                tri[k + 1][k + 2] = rest[1][1]
            else:
                fir = np.zeros((2, 2))
                fir = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).rd_to_b(fir)
                fir[0][0] = c
                fir[0][1] = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).simple_inverse_signal(s)
                fir[1][0] = s
                fir[1][1] = c
                sec = np.zeros((2, 1))
                sec = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).rd_to_b(sec)
                sec[0][0] = tri[n - 2][n - 1]
                sec[1][0] = tri[n - 1][n - 1]
                rest = self.simple_dot(fir, sec,c_operate)
                tri[n - 2][n - 1] = rest[0][0]
                tri[n - 1][n - 1] = rest[1][0]
        return tri, U, V

    def QR_zero_shift_once_iteration(self, tri, U, V, i_min, i_max,c_operate): # 已替换
        n = i_max + 1
        oldc = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(1)
        x = tri[i_min][i_min]
        y = tri[i_min][i_min + 1]  #########
        for i in range(i_min, n - 1, 1):
            c, s, r = self.csr(x, y,c_operate)
            V = self.updatecsvv(c, s, V, i,c_operate)
            if (i != i_min):
                tri[i - 1][i] = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).simple_inverse_signal(self.simple_time(olds, r,c_operate))
            x = self.simple_time(oldc, r,c_operate)
            y = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).simple_inverse_signal(self.simple_time(tri[i + 1][i + 1], s,c_operate))
            h = self.simple_time(tri[i + 1][i + 1], c,c_operate)
            c, s, r = self.csr(x, y,c_operate)
            U = self.updatecsvv(c, s, U, i,c_operate)
            # print('uuuu', U)
            tri[i][i] = r
            x = h
            if (i != n - 2):
                y = tri[i + 1][i + 2]
            oldc = c
            olds = s
        tri[n - 2][n - 1] = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).simple_inverse_signal(self.simple_time(h, s,c_operate))
        tri[n - 1][n - 1] = self.simple_time(h, c,c_operate)
        return tri, U, V

    def Bidiagonal_Singular_Value_Decomposition(self,Bi, U, V,c_operate):# 已替换
        epcl = 2 ** (-(c_operate[1]-1 ) * self.f_bits)
        underflow = 2 ** (-(c_operate[1]-1) * self.f_bits)  # 向下越界限，即机器可识别的最小正数
        tol = 100 * epcl
        it_num = 0
        old_i_low = -1
        old_i_high = -1
        tr=decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).rb_to_d(Bi)
        (m,n)=np.array(tr).shape
        fudge = min(m, n)
        maxit = 3 * n ** 2
        low_bound_sigma, lad, u = self.estimation_lower_bound_minimum_singular_value1(Bi,c_operate)
        low_bound_sigma=decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(low_bound_sigma)

        max_bound_sigma = abs(np.array(tr)).max()
        thresh = max(tol * low_bound_sigma, maxit * underflow)
        while (1):
            it_num += 1
            if it_num > 2000:
                raise EOFError('too much iteration!!')
            # print('it_num',it_num)
            i_upper_limit = n - 1
            # print('num', it_num)

            while (decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(self.simple_abs((Bi[i_upper_limit - 1][ i_upper_limit]))) <= thresh and i_upper_limit >= 1):
                # print('thresh',thresh)
                # print('now',self.b_to_d(self.simple_abs((Bi[i_upper_limit - 1][ i_upper_limit]))))
                i_upper_limit -= 1

            if i_upper_limit == 0:
                break
            i_low = i_upper_limit - 1
            if i_low!=0:
                while (decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(self.simple_abs((Bi[i_low - 1][ i_low]))) > thresh and i_upper_limit >= 1):
                    i_low -= 1
                    if i_low==0:
                        break
            # print('min,max', i_low, i_upper_limit)
            if i_upper_limit == i_low + 1:
                B_temp, U_temp1, U_temp2, V= self.mul_22_submatrix(
                    self.sector(Bi, [i_low, i_upper_limit], [i_low, i_upper_limit]),
                    self.sector(U, [0, len(U) - 1], [i_low, i_low]),
                    self.sector(U, [0, len(U) - 1], [i_upper_limit, i_upper_limit]),V,i_low,c_operate)
                Bi[i_low][i_low]=B_temp[0][0]
                Bi[i_low][i_low+1] = B_temp[0][1]
                Bi[i_low+1][i_low] = B_temp[1][0]
                Bi[i_low+1][i_low+1] = B_temp[1][1]
                for k in range(len(U)):
                    U[k][i_low]=U_temp1[k][0]
                    U[k][i_upper_limit]=U_temp2[k][0]
                break #??这里有问题？？？
            if (old_i_low != i_low or old_i_high != i_upper_limit):
                if abs(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(Bi[i_low][ i_low])) >= abs(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(Bi[i_upper_limit][ i_upper_limit])):  # 此处还要加判定条件：此次迭代的矩阵和上次的不同才需要重新设置方向
                    direction = 'down'
                else:
                    direction = 'up'
            old_i_low = i_low
            old_i_high = i_upper_limit
            direction = 'down'
            # print('direction', direction)
            # Apply convergence criteria__还没写完
            if direction == 'down':
                if abs(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(self.simple_divide(Bi[i_upper_limit - 1][ i_upper_limit] , lad[i_upper_limit],c_operate))) <= tol:
                    Bi[i_upper_limit - 1][ i_upper_limit] = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)
                for k in range(i_upper_limit-1):
                    if decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(self.simple_abs(self.simple_divide(Bi[k][k+1],u[k],c_operate)))<=tol:###这里可能有点问题
                        Bi[k][k + 1]=decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)
            # Compute shift
            # print('fudge*tol*low_bound_sigma/max_bound_sigma', fudge * tol * low_bound_sigma / max_bound_sigma)
            if fudge * tol * low_bound_sigma / max_bound_sigma <= epcl:
                shift = 0
            else:
                if direction == 'down':
                    s = Bi[i_upper_limit][ i_upper_limit]
                    d = self.simple_minus([self.simple_add(
                        [self.simple_power_two(Bi[i_upper_limit - 1][i_upper_limit - 1],c_operate), self.simple_power_two(Bi[
                                                                                                                    i_upper_limit - 2][
                                                                                                                    i_upper_limit - 1],c_operate)],c_operate),
                                           self.simple_add([self.simple_power_two(Bi[i_upper_limit][i_upper_limit],c_operate),
                                                            self.simple_power_two(
                                                                Bi[i_upper_limit - 1][i_upper_limit],c_operate)],c_operate)],c_operate)
                    d=self.simple_divide(d,decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(2),c_operate)
                    shift = self.simple_minus([self.simple_add(
                        [self.simple_power_two(Bi[i_upper_limit][i_upper_limit],c_operate), self.simple_power_two(Bi[
                                                                                                            i_upper_limit - 1][
                                                                                                            i_upper_limit],c_operate),
                         d],c_operate), self.simple_time(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(np.sign(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(d))), self.simple_sqrt(self.simple_add(
                        [self.simple_power_two(d,c_operate),
                         self.simple_time(self.simple_power_two(Bi[i_upper_limit - 1][i_upper_limit - 1],c_operate),
                                          self.simple_power_two(Bi[i_upper_limit - 1][i_upper_limit],c_operate),c_operate)],c_operate),c_operate),c_operate)],c_operate)
                else:
                    s = Bi[i_low, i_low]
                    d = ((Bi[i_low + 1, i_low + 1] ** 2 + Bi[i_low + 1, i_low + 2] ** 2) - (
                            Bi[i_low, i_low] ** 2 + Bi[i_low, i_low + 1] ** 2)) / 2
                    shift = (Bi[i_low, i_low] ** 2 + Bi[i_low, i_low + 1] ** 2) + d - np.sign(d) * np.sqrt(
                        d ** 2 + Bi[i_low + 1, i_low + 1] ** 2 * Bi[i_low, i_low + 1] ** 2)
                if decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(self.simple_power_two(self.simple_divide(shift, s,c_operate) ,c_operate)) <= epcl:
                    shift = 0
            #print('shift', self.b_to_d(shift))
            # / *Perform QR iteration * /
            if shift != 0:
                if direction == 'down':
                    Bi, U ,V= self.QR_Wilkinson_shift_Iteration_once(Bi, U, V, i_low, i_upper_limit, shift,c_operate)
                    if decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(self.simple_abs(Bi[i_upper_limit-1][i_upper_limit]))<=thresh:
                        Bi[i_upper_limit-1][i_upper_limit]=decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)
                else:
                    Bi, U, V = self.QR_Wilkinson_shift_Iteration_once(Bi, U, V, i_low, i_upper_limit, shift,c_operate)
            else:
                if direction == 'down':
                    Bi, U,V= self.QR_zero_shift_once_iteration(Bi, U, V,  i_low, i_upper_limit,c_operate)
                    if decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(self.simple_abs(Bi[i_upper_limit-1][i_upper_limit]))<=thresh:
                        Bi[i_upper_limit-1][i_upper_limit]=decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).d_to_b(0)
                else:
                    Bi, U = QR_zero_shift_once_iteratin_upward(Bi, U, i_low, i_upper_limit)
        Vt=self.conjt(V)
        return Bi, U, Vt

    def my_SVD2(self, mat,c_operate, test_mode=True): # 替换
        #print("my_SVD2__________________")
        #print("cl_House",cl_House)

        U, A, V = self.house_bidiag_explicit_UV(mat,c_operate)
        A = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).complex_to_real(A)
        A, U, Vt = self.Bidiagonal_Singular_Value_Decomposition(A, U, V,c_operate)
        if test_mode == True:
            for i in range(len(A[0])):
                if decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(A[i][i]).real < 0:
                    A[i][i] = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).simple_inverse_signal(A[i][i])
                    for k in range(len(U)):
                        U[k][i] = decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).inverse_signal(U[k][i])
            for i in range(len(A[0]) - 1):
                for j in range(i + 1, len(A[0])):
                    if (abs(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(A[i][i])) < abs(decimal_HPIEEE754_transform(c_operate[1], self.e_bits, self.f_bits).b_to_d(A[j][j]))):
                        tmp1 = A[i][i]
                        A[i][i] = A[j][j]
                        A[j][j] = tmp1
                        for k in range(len(U)):
                            tmp2 = U[k][i]
                            U[k][i] = U[k][j]
                            U[k][j] = tmp2
                        for k in range(len(Vt[0])):
                            tmp2 = Vt[i][k]
                            Vt[i][k] = Vt[j][k]
                            Vt[j][k] = tmp2
        return U, A, Vt

    def HSVD_L(self, mat,c_operate):   #已替换
        # ff, gg = mat.shape
        ff = len(mat)
        gg = len(mat[0])
        if ff >= gg:
            U2, S2, V2 = self.my_SVD2(mat, c_operate,True)
        else:
            V2, S2, U2 = self.my_SVD2(self.T(mat),c_operate,True)
            V2 = self.T(V2)
            U2 = self.T(U2)
        return U2, S2, V2

def test(calculator, num_test, c_operate, matrixsize, num_samples):
    relative_error_list = []
    time_consume = 0
    print("num_samples",num_samples)
    if calculator == 'simple add':
        for i in range(num_samples):
            print("num_samples", i)
            d1 = random.uniform(0, 10) ** random.uniform(-20, 20)
            d2 = random.uniform(0, 10) ** random.uniform(-20, 20)
            time_set = time.time()
            BFP_res = bnpc.simple_add([d_H_trans[c_operate[0]].d_to_b(d1), d_H_trans[c_operate[0]].d_to_b(d2)],c_operate)
            time_consume += time.time() - time_set
            d1 = decimal.Decimal(d1)
            d2 = decimal.Decimal(d2)
            accurate_res = d1 + d2
            relative_error = (decimal.Decimal(str(d_H_trans[c_operate[0]].b_to_d(BFP_res))) - accurate_res) / accurate_res
            relative_error_list.append(relative_error)

    if calculator == 'simple minus':
        for i in range(num_samples):
            d1 = random.uniform(0, 10) ** random.uniform(-20, 20)
            d2 = random.uniform(0, 10) ** random.uniform(-20, 20)
            time_set = time.time()
            BFP_res = bnpc.simple_minus([d_H_trans[c_operate[0]].d_to_b(d1), d_H_trans[c_operate[0]].d_to_b(d2)],c_operate)
            time_consume += time.time() - time_set
            d1 = decimal.Decimal(d1)
            d2 = decimal.Decimal(d2)
            accurate_res = d1 - d2
            relative_error = (decimal.Decimal(str(d_H_trans[c_operate[0]].b_to_d(BFP_res))) - accurate_res) / accurate_res
            relative_error_list.append(relative_error)



    if calculator == 'simple time':
        for i in range(num_samples):
            d1 = random.uniform(0, 10) ** random.uniform(-20, 20)
            d2 = random.uniform(0, 10) ** random.uniform(-20, 20)
            time_set = time.time()
            BFP_res = bnpc.simple_time(d_H_trans[c_operate[0]].d_to_b(d1), d_H_trans[c_operate[0]].d_to_b(d2),c_operate)
            time_consume += time.time() - time_set
            d1 = decimal.Decimal(d1)
            d2 = decimal.Decimal(d2)
            accurate_res = d1 * d2
            relative_error = (decimal.Decimal(str(d_H_trans[c_operate[0]].b_to_d(BFP_res))) - accurate_res) / accurate_res
            relative_error_list.append(relative_error)

    if calculator == 'simple divide':
        for i in range(num_samples):
            d1 = random.uniform(0, 10) ** random.uniform(-20, 20)
            d2 = random.uniform(0, 10) ** random.uniform(-20, 20)
            time_set = time.time()
            BFP_res = bnpc.simple_divide(d_H_trans[c_operate[0]].d_to_b(d1), d_H_trans[c_operate[0]].d_to_b(d2),c_operate)
            time_consume += time.time() - time_set
            d1 = decimal.Decimal(d1)
            d2 = decimal.Decimal(d2)
            accurate_res = d1 / d2
            relative_error = (decimal.Decimal(str(d_H_trans[c_operate[0]].b_to_d(BFP_res))) - accurate_res) / accurate_res
            relative_error_list.append(relative_error)
            # if relative_error > 10 ** -5:
            #     print('d1 ', decimal.Decimal(d1))
            #     print('d2 ', decimal.Decimal(d2))
            #     print('accurate_res ', decimal.Decimal(accurate_res))
            #     print("decimal.Decimal(str(d_H_trans[c_operate[0]].b_to_d(BFP_res)))",d_H_trans[c_operate[0]].b_to_d(BFP_res))
            #     raise EOFError('accuracy not meet!')

    if calculator == 'simple sqrt':
        for i in range(num_samples):
            d1 = random.uniform(0, 10) ** random.uniform(-20, 20)
            time_set = time.time()
            BFP_res = bnpc.simple_sqrt(d_H_trans[c_operate[0]].d_to_b(d1),c_operate)
            time_consume += time.time() - time_set
            d1 = decimal.Decimal(d1)
            accurate_res = math.sqrt(d1)
            accurate_res=decimal.Decimal(accurate_res)
            relative_error = (decimal.Decimal(str(d_H_trans[c_operate[0]].b_to_d(BFP_res))) - accurate_res) / accurate_res
            relative_error_list.append(relative_error)

    if calculator == 'complex add':

        for i in range(num_samples):

            d1 = random.uniform(0, 10) ** random.uniform(-20, 20) + 1j*random.uniform(0, 10) ** random.uniform(-20, 20)
            d2 = random.uniform(0, 10) ** random.uniform(-20, 20) + 1j*random.uniform(0, 10) ** random.uniform(-20, 20)
            time_set = time.time()
            BFP_res = bnpc.add([d_H_trans[c_operate[0]].complex_d_to_b(d1), d_H_trans[c_operate[0]].complex_d_to_b(d2)],c_operate)
            time_consume += time.time() - time_set
            accurate_res = d1 + d2
            relative_error = abs(((d_H_trans[c_operate[0]].complex_b_to_d(BFP_res))) - accurate_res) / abs(accurate_res)
            print("relative_error", relative_error)
            relative_error_list.append(relative_error)

    if calculator == 'complex minus':
        for i in range(num_samples):
            d1 = random.uniform(0, 10) ** random.uniform(-20, 20) + 1j*random.uniform(0, 10) ** random.uniform(-20, 20)
            d2 = random.uniform(0, 10) ** random.uniform(-20, 20) + 1j*random.uniform(0, 10) ** random.uniform(-20, 20)
            time_set = time.time()
            BFP_res = bnpc.minus(d_H_trans[c_operate[0]].complex_d_to_b(d1), d_H_trans[c_operate[0]].complex_d_to_b(d2),c_operate)
            time_consume += time.time() - time_set
            accurate_res = d1 - d2
            relative_error = abs(((d_H_trans[c_operate[0]].complex_b_to_d(BFP_res))) - accurate_res) / abs(accurate_res)
            relative_error_list.append(relative_error)



    if calculator == 'complex time':
        for i in range(num_samples):
            d1 = random.uniform(0, 10) ** random.uniform(-20, 20) + 1j*random.uniform(0, 10) ** random.uniform(-20, 20)
            d2 = random.uniform(0, 10) ** random.uniform(-20, 20) + 1j*random.uniform(0, 10) ** random.uniform(-20, 20)
            time_set = time.time()
            BFP_res = bnpc.time(d_H_trans[c_operate[0]].complex_d_to_b(d1), d_H_trans[c_operate[0]].complex_d_to_b(d2),c_operate)
            time_consume += time.time() - time_set
            accurate_res = d1 * d2
            relative_error = abs(((d_H_trans[c_operate[0]].complex_b_to_d(BFP_res))) - accurate_res) / abs(accurate_res)
            relative_error_list.append(relative_error)

    if calculator == 'complex divide':
        for i in range(num_samples):
            d1 = random.uniform(0, 10) ** random.uniform(-20, 20) + 1j*random.uniform(0, 10) ** random.uniform(-20, 20)
            d2 = random.uniform(0, 10) ** random.uniform(-20, 20) + 1j*random.uniform(0, 10) ** random.uniform(-20, 20)
            time_set = time.time()
            BFP_res = bnpc.divide(d_H_trans[c_operate[0]].complex_d_to_b(d1), d_H_trans[c_operate[0]].complex_d_to_b(d2),c_operate)
            time_consume += time.time() - time_set
            accurate_res = d1 / d2
            relative_error = abs(((d_H_trans[c_operate[0]].complex_b_to_d(BFP_res))) - accurate_res) / abs(accurate_res)
            relative_error_list.append(relative_error)
            # if relative_error > 10 ** -5:
            #     print('d1 ', decimal.Decimal(d1))
            #     print('d2 ', decimal.Decimal(d2))
            #     print('accurate_res ', decimal.Decimal(accurate_res))
            #     print("decimal.Decimal(str(d_H_trans[c_operate[0]].b_to_d(BFP_res)))",d_H_trans[c_operate[0]].b_to_d(BFP_res))
            #     raise EOFError('accuracy not meet!')

    if calculator == 'complex sqrt':
        for i in range(num_samples):
            d1 = random.uniform(0, 10) ** random.uniform(-20, 20) + 1j * random.uniform(0, 10) ** random.uniform(-20,20)
            time_set = time.time()
            print(" d1", d1)
            BFP_res = bnpc.sqrt(d_H_trans[c_operate[0]].complex_d_to_b(d1),c_operate)
            time_consume += time.time() - time_set
            accurate_res = cmath.sqrt(complex(d1.real,d1.imag))
            #print(cmath.sqrt(complex(d1.real,d1.imag)))
            relative_error = abs(((d_H_trans[c_operate[0]].complex_b_to_d(BFP_res))) - accurate_res) / abs(accurate_res)
            relative_error_list.append(relative_error)
    if calculator == 'matrix add':
        for i in range(num_samples):
            print("num_samples", i)
            A = np.random.randn(matrixsize, matrixsize) + 1j * np.random.randn(matrixsize, matrixsize)
            B = np.random.randn(matrixsize, matrixsize) + 1j * np.random.randn(matrixsize,matrixsize)
            A_BFP = d_H_trans[c_operate[0]].cd_to_b(A)
            B_BFP = d_H_trans[c_operate[0]].cd_to_b(B)
            add_standard = A+B
            time_set = time.time()
            add_BFP = bnpc.A_add_B(A_BFP, B_BFP, c_operate)
            time_consume += time.time() - time_set
            add_test = np.array(d_H_trans[c_operate[1]].cb_to_d(add_BFP))
            # 测试
            accurate_res = add_test - add_standard
            relative_error = np.linalg.norm(accurate_res)/np.linalg.norm(add_standard)
            relative_error_list.append(relative_error)
    if calculator == 'matrix minus':
        for i in range(num_samples):
            print("num_samples", i)
            A = np.random.randn(matrixsize, matrixsize) + 1j * np.random.randn(matrixsize, matrixsize)
            B = np.random.randn(matrixsize, matrixsize) + 1j * np.random.randn(matrixsize,matrixsize)
            minus_standard = A-B
            A_BFP = d_H_trans[c_operate[0]].cd_to_b(A)
            B_BFP = d_H_trans[c_operate[0]].cd_to_b(B)
            time_set = time.time()
            minus_BFP = bnpc.A_minus_B(A_BFP, B_BFP, c_operate)
            time_consume += time.time() - time_set
            minus_test = np.array(d_H_trans[c_operate[1]].cb_to_d(minus_BFP))
            # 测试
            accurate_res = minus_test - minus_standard
            relative_error = np.linalg.norm(accurate_res)/np.linalg.norm(minus_standard)
            relative_error_list.append(relative_error)
    if calculator == 'matrix time':
        for i in range(num_samples):
            print("num_samples", i)
            A = np.random.randn(matrixsize, matrixsize) + 1j * np.random.randn(matrixsize, matrixsize)
            B = np.random.randn(matrixsize, matrixsize) + 1j * np.random.randn(matrixsize,matrixsize)
            time_standard = np.dot(A, B)
            A_BFP = d_H_trans[c_operate[0]].cd_to_b(A)
            B_BFP = d_H_trans[c_operate[0]].cd_to_b(B)
            time_set = time.time()
            time_BFP = bnpc.dot(A_BFP, B_BFP, c_operate)
            time_consume += time.time() - time_set
            time_test = np.array(d_H_trans[c_operate[1]].cb_to_d(time_BFP))
            accurate_res = time_test - time_standard
            relative_error = np.linalg.norm(accurate_res)/np.linalg.norm(time_standard)
            relative_error_list.append(relative_error)
    if calculator == 'matrix svd':
        for i in range(num_samples):
            print('num_samples',i)
            A = np.random.randn(matrixsize, matrixsize) + 1j * np.random.randn(matrixsize, matrixsize)
            U1, S1, V1 = np.linalg.svd(A, 1, 1)
            S_standard = np.diag(S1)
            A_BFP = d_H_trans[c_operate[1]].cd_to_b(A)
            time_set = time.time()
            U_BFP, S_BFP, V_BFP = bnpc.HSVD_L(A_BFP, c_operate)
            time_consume += time.time() - time_set
            S_test =np.array(d_H_trans[c_operate[1]].rb_to_d(S_BFP))
            accurate_res = S_test - S_standard
            relative_error = np.linalg.norm(accurate_res)/np.linalg.norm(S_standard)
            relative_error_list.append(relative_error)
    if calculator == 'matrix inv':
        for i in range(num_samples):
            print("num_samples", i)
            A = np.random.randn(matrixsize, matrixsize) + 1j * np.random.randn(matrixsize, matrixsize)
            inv_standard = np.linalg.inv(A)
            A_BFP = d_H_trans[c_operate[0]].cd_to_b(A)
            time_set = time.time()
            inv_BFP= bnpc.complex_inv(A_BFP, c_operate)
            time_consume += time.time() - time_set
            inv_test =np.array(d_H_trans[c_operate[1]].cb_to_d(inv_BFP))
            accurate_res = inv_test - inv_standard
            relative_error = np.linalg.norm(accurate_res)/np.linalg.norm(inv_standard)
            relative_error_list.append(relative_error)
            if relative_error > 10 ** (-3):
                print("relative_error",relative_error)
                print("A_BFP",A_BFP)

    # print('complete!')
    # print('max relative error: ', np.max(relative_error_list))
    # print('min relative error: ', np.min(relative_error_list))
    # print('avg relative error: ', np.mean(relative_error_list))
    # print('var relative error: ', np.var(relative_error_list))
    # print('time consume: ', time_consume)
    time_consume = time_consume/num_samples
    return [np.max(relative_error_list), np.min(relative_error_list),
            abs(np.mean(relative_error_list)), np.var(relative_error_list), time_consume]


#
#
# def monte_test2(test_list, num_multiplies, num_samples,blocknum, matrixsize):#测试算子叠加误差
#     error_list = []
#     for i in range(len(test_list)):
#         error_list.append([])
#         for j in range(num_multiplies):
#             error_list[-1].append([])
#     for k in range(num_samples):
#         test_list_decimal = []
#         test_list_BFP = []
#         for i in range(num_multiplies):
#             flag = random.random()
#             d = np.random.randn(matrixsize, matrixsize) + 1j * np.random.randn(matrixsize, matrixsize)  # 生成测试矩阵
#             d_BFP = d_H_trans[blocknum].cd_to_b(d)
#             test_list_BFP.append(d_BFP)
#             #print('b: ', (test_list_BFP[-1]))
#             test_list_decimal.append(d)
#         BFP_res = []
#         decimal_res = []
#         BFP_store = []
#         for i in range(len(test_list)):
#             BFP_res.append(test_list_BFP[0])
#             BFP_store.append(test_list_BFP[0])
#             decimal_res.append(test_list_decimal[0])
#         for i in range(len(test_list)):
#             error0 = d_H_trans[blocknum].cb_to_d(BFP_res[0]) - decimal_res[0]
#             error_list[i][0].append(np.linalg.norm(error0) / np.linalg.norm(decimal_res[0]))
#         for i in range(num_multiplies - 1):
#             #print(test_list_decimal)
#             count = 0
#
#             if 'matrix add' in test_list:
#                 print("count", count)
#                 decimal_res[count] = decimal_res[count] + test_list_decimal[i + 1]
#                 BFP_res[count] = bnpc.A_add_B(BFP_res[count], test_list_BFP[i + 1], blocknum)
#                 error0 = d_H_trans[blocknum].cb_to_d(BFP_res[count]) - decimal_res[count]
#                 error_list[count][i + 1].append(np.linalg.norm(error0) / np.linalg.norm(decimal_res[count]))
#                 count += 1
#             if 'matrix minus' in test_list:
#                 print("count", count)
#                 decimal_res[count] = decimal_res[count] + test_list_decimal[i + 1]
#                 BFP_res[count] = bnpc.A_minus_B(BFP_res[count], test_list_BFP[i + 1], blocknum)
#                 error0 = d_H_trans[blocknum].cb_to_d(BFP_res[count]) - decimal_res[count]
#                 error_list[count][i + 1].append(np.linalg.norm(error0) / np.linalg.norm(decimal_res[count]))
#                 count += 1
#             if 'matrix time' in test_list:
#                 print("count", count)
#                 decimal_res[count] = np.dot(decimal_res[count], test_list_decimal[i + 1])
#                 BFP_res[count] = bnpc.dot(BFP_res[count], test_list_BFP[i + 1], blocknum)
#                 error0 = d_H_trans[blocknum].cb_to_d(BFP_res[count]) - decimal_res[count]
#                 error_list[count][i + 1].append(np.linalg.norm(error0) / np.linalg.norm(decimal_res[count]))
#                 count += 1
#
#             if 'simple time half' in test_list:
#                 decimal_res[count] = decimal_res[count] * test_list_decimal[i + 1]
#                 BFP_res[count] = bnp.simple_time_half(BFP_res[count], test_list_BFP[i + 1])
#                 error_list[count][i + 1].append((decimal.Decimal(str(d_H_trans.b_to_d(BFP_res[count])))
#                                                  - decimal_res[count]) / decimal_res[count])
#                 count += 1
#             if 'simple time all' in test_list:
#                 decimal_res[count] = decimal_res[count] * test_list_decimal[i + 1]
#                 BFP_res[count] = bnp.simple_time_all(BFP_res[count], test_list_BFP[i + 1])
#                 error_list[count][i + 1].append((decimal.Decimal(str(d_H_trans.b_to_d(BFP_res[count])))
#                                                  - decimal_res[count]) / decimal_res[count])
#                 count += 1
#             if 'time store' in test_list:
#                 #测试乘法对量化影响的误差
#                 decimal_res[count] = decimal_res[count] * test_list_decimal[i + 1]
#                 b1 = d_H_trans.b_to_d(test_list_BFP[i + 1])
#                 b2 = d_H_trans.b_to_d(BFP_store[count])
#                 b_time = b1 * b2
#                 BFP_store[count] = d_H_trans.d_to_b(b_time)
#                 error_list[count][i + 1].append((decimal.Decimal(str(d_H_trans.b_to_d(BFP_store[count])))
#                                                  - decimal_res[count]) / decimal_res[count])
#                 count += 1
#             if 'divide store' in test_list:
#                 # 测试量化误差对除法的影响
#                 decimal_b2 = decimal_res[count]
#                 decimal_res[count] = decimal_res[count] / test_list_decimal[i + 1]
#                 b1 = d_H_trans.b_to_d(test_list_BFP[i + 1])          #b1 = d_H_trans.b_to_d(test_list_BFP[i + 1])
#                 print('b1', b1)
#                 b2 = d_H_trans.b_to_d(BFP_store[count])
#                 b_divide = b2 / b1
#                 BFP_store[count] = d_H_trans.d_to_b(b_divide)
#                 error_list[count][i + 1].append((decimal.Decimal(str(d_H_trans.b_to_d(BFP_store[count])))
#                                                  - decimal_res[count]) / decimal_res[count])
#                 count += 1
#             if 'simple divide' in test_list:
#                 decimal_res[count] = decimal_res[count] / test_list_decimal[i + 1]
#                 BFP_res[count] = bnp.simple_divide(BFP_res[count], test_list_BFP[i + 1])
#                 error_list[count][i + 1].append((decimal.Decimal(str(d_H_trans.b_to_d(BFP_res[count])))
#                                                  - decimal_res[count]) / decimal_res[count])
#                 count += 1
#
#     mean_list = []
#     var_list = []
#     for i in range(len(test_list)):
#         mean_list.append([])
#         var_list.append([])
#         for j in range(num_multiplies):
#             tmp = np.array(error_list[i][j])
#             mean_list[i].append(np.mean(tmp))
#             var_list[i].append(np.var(tmp))
#     return mean_list, var_list
#
#
def relation_bytes_accuracy(calculator, num_test, matrixsize, start, end, step, test_type, fixed_cl,num_samples):
    statistics = [[], [], [], []]
    c_operate = copy.copy(fixed_cl)
    for var_cl in range(start, end, step):
        print("----------", var_cl, '------', calculator)
        if test_type == 'block of time':
            c_operate[3] = var_cl
            tmp = test(calculator, num_test, c_operate, matrixsize,num_samples)
            print("len(tmp)",len(tmp))
            for i in range(len(tmp)):
                statistics[i].append(tmp[i])

        if test_type == 'block of add and minus':
            c_operate[1] = var_cl
            c_operate[2] = var_cl
            tmp = test(calculator, num_test, c_operate, matrixsize,num_samples)
            print("len(tmp)",len(tmp))
            for i in range(len(tmp)):
                statistics[i].append(tmp[i])

        if test_type == 'block of divide':
            c_operate[4] = var_cl
            tmp = test(calculator, num_test, c_operate, matrixsize,num_samples)
            print("len(tmp)",len(tmp))
            for i in range(len(tmp)):
                statistics[i].append(tmp[i])

        if test_type == 'block of sqrt':
            c_operate[5] = var_cl
            tmp = test(calculator, num_test, c_operate, matrixsize,num_samples)
            print("len(tmp)",len(tmp))
            for i in range(len(tmp)):
                statistics[i].append(tmp[i])
    return statistics

def relation_invariable_bytes_accuracy(calculator, num_test, matrixsize, start, end, step, test_type, fixed_cl,num_samples):
    statistics = [[], [], [], [],[]]
    for var_cl in range(start, end, step):
        print("----------", var_cl, '------', test_type)
        if test_type == 'all':
            c_operate = copy.copy(fixed_cl)
            c_operate[1] = var_cl
            c_operate[2] = var_cl
            c_operate[3] = var_cl
            c_operate[4] = var_cl
            c_operate[5] = var_cl
            c_operate[0] = var_cl
            print("fixed_cl", fixed_cl)
            print("c_operate", c_operate)
            tmp = test(calculator, num_test, c_operate, matrixsize,num_samples)
            print("len(tmp)",len(tmp))
            for i in range(len(tmp)):
                statistics[i].append(tmp[i])
    print("statistics",statistics)
    return statistics

def relation_operate_accuracy(calculator, num_test, matrixsize, start, end, step, test_type, fixed_cl, num_samples):
    statistics = [[], [], [], []]
    for var_cl in range(start, end, step):
        print("----------", var_cl, '------', test_type)
        if test_type == 'block of time':
            c_operate = copy.copy(fixed_cl)
            c_operate[3] = var_cl
            print("fixed_cl", fixed_cl)
            print("c_operate", c_operate)
            tmp = test(calculator, num_test, c_operate, matrixsize,num_samples)

            print("len(tmp)",len(tmp))
            for i in range(len(tmp)):
                statistics[i].append(tmp[i])

        if test_type == 'block of add and minus':

            c_operate = copy.copy(fixed_cl)
            c_operate[1] = var_cl
            c_operate[2] = var_cl
            print("fixed_cl", fixed_cl)
            print("c_operate", c_operate)
            tmp = test(calculator, num_test, c_operate, matrixsize,num_samples)
            print("len(tmp)",len(tmp))
            for i in range(len(tmp)):
                statistics[i].append(tmp[i])

        if test_type == 'block of divide':
            c_operate = copy.copy(fixed_cl)
            c_operate[4] = var_cl

            print("c_operate", c_operate)
            tmp = test(calculator, num_test, c_operate, matrixsize,num_samples)
            print("len(tmp)",len(tmp))
            for i in range(len(tmp)):
                statistics[i].append(tmp[i])

        if test_type == 'block of sqrt':
            c_operate = copy.copy(fixed_cl)
            c_operate[5] = var_cl
            print("fixed_cl", fixed_cl)
            print("c_operate", c_operate)
            tmp = test(calculator, num_test, c_operate, matrixsize,num_samples)
            print("len(tmp)",len(tmp))
            for i in range(len(tmp)):
                statistics[i].append(tmp[i])
        if test_type == 'all':
            c_operate = copy.copy(fixed_cl)
            c_operate[1] = var_cl
            c_operate[2] = var_cl
            c_operate[3] = var_cl
            c_operate[4] = var_cl
            c_operate[5] = var_cl
            c_operate[0] = var_cl
            print("fixed_cl", fixed_cl)
            print("c_operate", c_operate)
            tmp = test(calculator, num_test, c_operate, matrixsize,num_samples)
            print("len(tmp)",len(tmp))
            for i in range(len(tmp)):
                statistics[i].append(tmp[i])
    print(test_type,"statistics",statistics)
    return statistics

def relation_matrixsize_accuracy(calculator, num_test,blocknum, msizestart, msizeend,fixed_cl,num_samples,step):
    statistics = [[], [], [], [],[]]
    for msize in range(msizestart, msizeend, step):
        print("----------",msize,'------',calculator)
        tmp = test(calculator, num_test,fixed_cl,msize,num_samples)
        print("len(tmp)", len(tmp))
        for i in range(len(tmp)):
            statistics[i].append(tmp[i])
    return statistics

def plt_graph(form, x, ylist, title, xlabel, ylabel, label_list,x_log_flag, y_log_flag):
    num = []
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if form == 'plot':
        for i in range(len(ylist)):
            plt.plot(x, ylist[i], label=label_list)
    if form == 'scatter':
        for i in range(len(ylist)):
            plt.scatter(x, ylist[i], label=label_list)
    if x_log_flag == 1:
        plt.xscale('log')
    if y_log_flag == 1:
        plt.yscale('log')
    my_x_ticks = x  # 原始数据有13个点，故此处为设置从0开始，间隔为1
    plt.xticks(my_x_ticks)
    plt.legend(labels=label_list)


if __name__ == '__main__':
    a=1
    bnpc = bnpc(13, 8)
    d_H_trans = generate_d_H_trans(256, 13, 8)  # (256,8,4)
    # add_list = [[[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1], [1], [1], [1], [1], [1], [1], [1], [1], [0], [0], [0], [1], [0], [0]]]
    # bnpc.simple_add_same_signal(add_list, 8)

    ######## 测试基本实数算子的误差
    #fixed_cl =[8, 8, 8, 8, 8,8]
    # statistic_list = []
    # #test_list = ['simple add', 'simple minus', 'simple time', 'simple divide', 'simple sqrt']
    # #test_list = ['complex add', 'complex minus', 'complex time', 'complex divide']
    # test_list = ['matrix add', 'matrix minus', 'matrix time','matrix inv','matrix svd']
    # test_type = 'all'
    # matrixsize = 3
    # num_test = 300
    # step = 1
    # start = 4
    # end = 13
    # num_samples = 100
    # fixed_cl = [8, 8, 8, 8, 8, 8]
    # for i in range(len(test_list)):
    #     statistic_list.append(relation_invariable_bytes_accuracy(test_list[i], num_test, matrixsize, start, end, step, test_type, fixed_cl,num_samples))
    #     statistic_list[-1] = np.array(statistic_list[-1])
    # statistic_list_rearrange = []
    # for i in range(len(statistic_list[0])):
    #     statistic_list_rearrange.append([])
    #     for j in range(len(statistic_list)):
    #         statistic_list_rearrange[-1].append(statistic_list[j][i])
    # plt.figure()
    # plt_graph('plot', list(range(start, end, step)), statistic_list_rearrange[2],
    #           'NMSE', 'number of blocks', 'E(error)', test_list, 0 , 1)
    # plt.show()
    # plt.figure()
    # plt_graph('plot', list(range(start, end, step)), statistic_list_rearrange[3],
    #           'Var','number of blocks', 'Var(error)', test_list, 0 , 1)
    # plt.show()
    # plt_graph('plot', list(range(start, end)), statistic_list_rearrange[4],
    #           'Time complexity of 1 times of operations', 'number of blocks', 'CPU time(s)', test_list, 0, 0)
    # plt.show()

    ####### 固定加法位数测试矩阵大小对误差的影响
    fixed_cl = [8, 16, 16, 16, 16, 16]
    statistic_list = []
    num_samples = 100
    #test_list = ['matrix add', 'matrix minus', 'matrix time','matrix inv','matrix svd']
    test_list = ['matrix inv']
    num_test = 5
    msizestart = 2
    step = 2
    msizeend = 10
    blocknum = 5
    for i in range(len(test_list)):
        statistic_list.append(relation_matrixsize_accuracy(test_list[i], num_test, blocknum, msizestart, msizeend,fixed_cl,num_samples,step))
        statistic_list[-1] = np.array(statistic_list[-1])
    statistic_list_rearrange = []
    for i in range(len(statistic_list[0])):
        statistic_list_rearrange.append([])
        for j in range(len(statistic_list)):
            statistic_list_rearrange[-1].append(statistic_list[j][i])
    plt.figure()
    plt_graph('plot', list(range(msizestart, msizeend, step)), statistic_list_rearrange[2],
              'NMSE', 'size of matrix', 'E(error)', test_list,0,1)
    plt.show()
    plt.figure()
    print("statistic_list_rearrange[3]",statistic_list_rearrange[3])
    plt_graph('plot', list(range(msizestart, msizeend, step)), statistic_list_rearrange[3],
              'Var', 'size of matrix', 'Var(error)', test_list,0,1)
    plt.show()

    ################################# 测试级联操作的误差
    # num_operations = 10
    # num_samples = 100
    # matrixsize =4
    # blocknum = 16
    # test_list = ['matrix add','matrix minus','matrix time']
    # mean_list, var_list = monte_test2(test_list, num_operations, num_samples,blocknum,matrixsize)
    # plt.figure()
    # plt_graph('plot', range(1, num_operations + 1), mean_list, 'chain mean',
    #           'time of operations', 'error', test_list)
    # plt.figure()
    # plt_graph('plot', range(1, num_operations + 1), var_list, 'chain var',
    #           'time of operations', 'error', test_list)
    #
    # plt.figure()
