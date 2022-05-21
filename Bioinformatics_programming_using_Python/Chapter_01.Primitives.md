## Introduction
to use the Python programming language to facilitate and automate the wide variety of data manipulation tasks encountered in life science research and development.
#### about bioinformatics
- *computaional biology*
- *Software development*
- *Life science research and development*

the book focuses on particular data management and manipulation tasks.

broadly speaking, programming languages embody combinations of four paradigms.
- procedural
- declarative
- functional
- object oriented

## primitives
**type** and **value**  *vs* **class** and **object**

### simple values
logical(Boolean), integer, float, and string

#### expressions
an *operator* is a symbol that indicates a calculation usng one or more operands. The combination of the operator and its operand is an expression
- numeriac operators:
- logical operations:      
  six comparison operators ==, !=, <, <=, >, >=
- string operations:
  four binary operators that act on strings: in, not in, + and * , and subscription (extracts character sunstring of a string), sclicing (extracts a series of characters from a string)
#### calls:  A call is a kind of expression.
*function calls* The function is called, does something, then returns a value.
```
len(arg)
print()
input(string)
abs(value)
max(args...)
min(args...)
str(arg)
int(arg)
float(arg)
bool(arg)
help()
```
*method calls* calling a method is just like calling a function, except taht the first argument goes before the function name, followed by a period.
```
string1.count(string2[,start[,end]])
string1.find(string2[,start[,end]])
string1.startswith(string2[,start[,end]])
string1.strip([string2])
string1.lstrip([string2])
string1.rstrip([string2])
```
#### compound expressions
