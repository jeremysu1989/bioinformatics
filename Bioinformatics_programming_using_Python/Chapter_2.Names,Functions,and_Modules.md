name values, define new functions, incorporate optional software from the python library

#### assigning names
a python program consists of a series of *statements*. statements do not produce values and cannot be used as expressions. ***An assignment statement binds a name to an object*
```
a = a + 1
a += 1
# argumented assignment statement
```

#### defining functions
New functions are defined with *function definition statements.* A definition is a compound statement, meaning it comprises more than one line of code. This first line of each cmpund statement ends with a colon.
```
def name(parameter-list):
    body
    return value
```
Do Nothing
```
def fn():
    pass
```

#### function parameters
```
def recognition_site(base_seq, recognition_seq):
    return base_seq.find(recognition_seq)
```


