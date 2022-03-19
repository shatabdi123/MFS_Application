# Python3 program to Convert a
# list to dictionary

def Convert(lst):
    res_dct = {lst[i]: 1 for i in range(0, len(lst), 1)}
    return res_dct


# Driver code
lst = ['a','b','c']
print(Convert(lst))