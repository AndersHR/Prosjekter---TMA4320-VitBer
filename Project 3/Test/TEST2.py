from numpy import array

def func(arr1, arr2):
    return -arr1**2/(arr1+arr2)

arr1, arr2 = array([1, 2 ,3]), array([5, 4, 7])
print (func(arr1, arr2))