def stringint(k, n):
    strint = str(k)
    res = '0' * (n - len(strint)) + strint
    return res
