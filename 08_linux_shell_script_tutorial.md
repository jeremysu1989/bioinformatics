## 1. Quick introduction to Linux
#### Unix philosophy
- do one thing and do it well
- everything is file
- small is beautiful
- store data and configuration in flat text files
- use shell scripts to incease leverage and portability
- chain programs together to complete complex task
- choose portability over efficiency
- keep it simple, stupid

#### What is a shell script or shell scripting
Each shell script consists of 
- shell keywords such asif..else, do..while
- shell commands such as pwd, test, echo
- linux binary commands such as who, free
- test processing utilities such as grep, awk,cut
- functions
- control flow statements such as if..then..else

Each script has purpose
- specific purpose
- act like a command
- script code usability

## 2. Getting started with shell programming
#### The role of shells in the linux enviroment
Shell is used for various prpose under lnux. Linux user enviroment is made of the following components.
- kernel - The core of linux operating syste,
- shell - provides an interface between the user and the kernel
- ternimal emulator - the xterm program is a terminal emulator for the X windows system
- linux desktop and window manager - linux desktop is collection of various software apps 

#### shebang
The #! syntas used in scripts to indicate an interpreter for execution under UNIX/linux operating systems
```
#!/bin/bash
OR
#!/usr/bin/perl
OR
#!/usr/bin/python
```

#### shell commnets
```
#!/bin/bash
# A Simple Shell Script To Get Linux Network Information
# Vivek Gite - 30/Aug/2009
echo "Current date : $(date) @ $(hostname)"
echo "Network configuration"
/sbin/ifconfig
```
multiple line comment
```
#!/bin/bash
echo "Adding new users to LDAP Server..."
<<COMMENT1
  Master LDAP server : dir1.nixcraft.net.in
  Add user to master and it will get sync to backup server too
  Profile and active directory hooks are below
COMMENT1
echo "Searching for user..."
```

## 3. The shell variables and environment
#### Variables in shell

#### Rules for naming variable name
*Do not put spaces on either side of the equal sign when assigning value to variable*
```
n=10
```

#### getting user input via keyboard
You can accept input from the keyboard and assign an input value to a user defined shell variable using read command
*read command syntxa*
```
read -p "prompt" variable1 varaiable2 variableN
```

#### path name expansion
```
{ pattern1, pattern2, patternN }
text{ pattern1, pattern2, patternN }
text1{ pattern1, pattern2, patternN }text2
command something/{ pattern1, pattern2, patternN }

echo file{1..5}.txt
```

#### changing bash prompt
```
# regular colors
local K="\[\033[0;30m\]" # black
local R="\[\033[0;31m\]" # red
local G="\[\033[0;32m\]" # green
local Y="\[\033[0;33m\]" # yellow
local B="\[\033[0;34m\]" # blue
local M="\[\033[0;35m\]" # magenta
local C="\[\033[0;36m\]" # cyan
local W="\[\033[0;37m\]" # white
```

## 4. Conditionals exections(Decision making)
#### bash strucutred language constructs
```
if today is friday
  execute tar command
otherwise
  print an error on screen.
```

#### if structures to execute code based on a condition
```
if conditon
then
  command1
  command2
  ...
  commandN
 fi
 ```
 OR
 ```
 if test var == value
 then
   command1
   command2
   ...
   commandN
 fi
 ```
#### conditional execution
```
command1 && command2
command1 || command2
```
#### logical not!
```
! expression
```
#### string comparison
```
string1 = string2
```
#### File attributes comparions
```
#!/bin/bash
FILEPATH="/Users/scarecrow/Temp/test"
DATAPATH="/Users/scarecrow/Data/seq_data"

read -p "Enter a temp folder name : "  DOMAIN

# Make sure we got the Input else die with an error on screen
[ -z $DOMAIN ] && { echo "Please enter a domain name. Try again!"; exit 1; }

# Alright, set some variable based upon $DOMAIN
OUT="$FILEPATH/$DOMAIN/today/logfolder"
TESTDATA="$DATAPATH/$DOMAIN/today"

# Die if configuration file exits...
[ -f $OUT ] && { echo "Configuration file '$OUT' exits for domain $FILEPATH."; exit 2; }

# Make sure configuration directory exists
[ ! -d $TESTDATA ] && mkdir -p $TESTDATA

mkdir -p $OUT
# Write a log file
touch $OUT/logfile.txt

result="$OUT/logfile.txt"

echo "$date" >> $result
echo "This is new." >> $result
echo "Newbee" >> $result
```

#### shell command ine parameters
- command line options
- options
- position paramenters
- flag
- swithes or swith
- command line arguments

#### How to use positional parameters
All command line parameters ( positional parameters ) are available via special shell variable $1, $2, $3,...,$9
```
#!/bin/bash
echo "The script name : $0"
echo "The value of the first argument to the script : $1"
echo "The value of the second argument to the script : $2"
echo "The value of the third argument to the script : $3"
echo "The number of arguments passed to the script : $#"
echo "The value of all command-line arguments (\$* version) : $*"
echo "The value of all command-line arguments (\$@ version) : $@"
```

#### The case statement
The case statement is good alternative to multilevel if-then-else-fi statement. It enable you to match several values against one variable. It is easier to read and write
```
case $variable-name in
    pattern1)
        command1
        ...
        ....
        commandN
        ;;
    pattern2)
        command1
        ...
        ....
        commandN
        ;;
    patternN)
        command1
        ...
        ....
        commandN
        ;;
    *)
esac
```


## 5. Bash Loops
#### the for loop statement
```
for var in item1 item2 ... itemN
do
        command1
        command2
        ....
        ...
        commandN
done
```
the for loop numerical explicit list syntax:
```
for var in list-of-values
do
        command1
        command2
        ....
        ...
        commandN
done
```
the foor loop explicit file list syntax:
```
for var in file1 file2 file3 fileN
do
        command1
        command2
        ....
        ...
        commandN
done
```
the for loop variable's content syntax:
```
for var in $fileNames
do
        command1
        command2
        ....
        ...
        commandN
done
```
the for loop command substitution syntax:
```
for var in $(Linux-command-name)
do
        command1
        command2
        ....
        ...
        commandN
done
```
the for loop explicit file list using hash arrry syntax:
```
# define an array
ArrayName=(~/.config/*.conf)
for var in "${ArrayName[@]}"
do
        command1 on $var
        command2
        ....
        ...
        commandN
done
```
the for loop three-expression syntax
```
for (( EXP1; EXP2; EXP3 ))
do
        command1
        command2
        command3
done
```
*examples*
```
#!/bin/bash
for i in 1 2 3 4 5
do
  echo "Welcome $i times."
done
```

```
#!/bin/bash
# A simple shell script to print list of cars
for car in bmw ford toyota nissan
do
  echo "Value of car is: $car"
done
```

```
#!/bin/bash
# A simple shell script to run commands
for command in date pwd df
do
  echo
  echo "*** The output of $command command >"
  #run command
  $command
  echo
done
```

```
#!/bin/bash
# A shell script to verify user password database
files="/etc/passwd /etc/group /etc/shadow /etc/gshdow"
for f in $files
do
  [ -f $f ] && echo "$f file found" || echo "*** Error - $f file missing."
done
```

#### nested for loop statement
```
#!/bin/bash
# A shell script to print each number five times.
for (( i = 1; i <= 5; i++ )) ### Outer for loop ###
do
    for (( j = 1 ; j <= 5; j++ )) ### Inner for loop ###
    do
        echo -n "$i "
    done
    echo "" #### print the new line ###
done
```
Chessboard example
```
#!/bin/bash
for (( i = 1; i <= 8; i++ )) ### Outer for loop ###
do
    for (( j = 1 ; j <= 8; j++ )) ### Inner for loop ###
    do
        total=$(( $i + $j)) # total
        tmp=$(( $total % 2)) # modulus
        # Find out odd and even number and change the color
        # alternating colors using odd and even number logic
        if [ $tmp -eq 0 ];
        then
            echo -e -n "\033[47m "
        else
            echo -e -n "\033[40m "
        fi
    done
    echo "" #### print the new line ###
done
```
#### The while loop statement
```
while [ condition ]
do
    command1
    command2
    ..
    ....
    commandN
done
```
To read a text file line-by-line, use the following syntax
```
while IFS= read -r line
do
    command1 on $line
    command2 on $line
    ..
    ....
    commandN
done < "/path/to/filename"
```
OR
```
while IFS= read -r field1 filed2 field3 ... fieldN
do
    command1 on $field1
    command2 on $field1 and $field3
    ..
    ....
    commandN on $field1 ... $fieldN
done < "/path/to dir/file name with space"
```
*examples*
```
#!/bin/bash
# set n to 1
n=1
# continue until $n equals 5
while [ $n -le 5 ]
do
    echo "Welcome $n times."
    n=$(( n+1 )) # increments $n
done
```
```
#!/bin/bash
n=1
while (( $n <= 5 ))
do
    echo "Welcome $n times."
    n=$(( n+1 ))
done
```
you can read a test file using read command and while loop as follows
```
#!/bin/bash
file=/etc/resolv.conf
while IFS= read -r line
do
    # echo line is stored in $line
    echo $line
done < "$file"

#outputs:
nameserver 127.0.0.1
nameserver 192.168.1.254
nameserver 4.2.2.1
```
you can store above output in two seperate fields as follows
```
#!/bin/bash
file=/etc/resolv.conf
while IFS= read -r f1 f2
do
    echo "field # 1 : $f1 ==> field #2 : $f2"
done < "$file"

#outputs
field # 1 : nameserver ==> field #2 : 127.0.0.1
field # 1 : nameserver ==> field #2 : 192.168.1.254
field # 1 : nameserver ==> field #2 : 4.2.2.1
```
##### command substitution
Command substitution is nothing but run a shell command and store it's output to a variable or display back using echo command. 
```
NOW=$(date)
echo "$NOW"
```
command substitution and shell loops
```
for f in $(ls /etc/*.conf)
do
    echo "$f"
done
```
## 6. Shell redirection
Changing the default path of input or output is called redirection
- In linux everything is a file
- Your hardware is also a file
  -   0 Input - Keyboard(stdin)
  -   1 Output - Screen(stdout)
  -   2 Error - Screen(stderr)

#### standard input
- standard input is the default input method
- it is denoted by zero number
- also known as stdin
- the default standard input is the keyboard
- < si input redirection symbol and syntax is 
```
command < filename
```

#### standard putput
- standard output is used by a command to writes its output
- the default is the screen
- it is denoted by one umber
- also known as stdout
- the default standard output is the screen
-  > is output redirection symbol and syntax is
```
command > output.file.name
```

#### standard error
- standard error is the default error output device, which is used to wirite all system error messages
- it is denoted by two number
- also know as stderr
- the default standard input is the screen or monitor
- 2> is input redirection symbol and syntax is
```
command 2> errors.txt
```

#### empty file creation
```
>newfile.name
```

### /dev/null discards unwantted output
All data written o a /dev/null or /dev/zero file is discarded by the system
```
command > /dev/null
```
### Here documents
```
command <<HERE
text1
text2
text3
$whatever
HERE
```
This type of redirection tells the shell to read input from the current source (HERE) until a line containg only word (HERE) is seen.


## 7. Pipes and filters
#### commonly used filter commands
```
awk.   cut.   grep.   gzip.   head.   paste.   perl    
sed.   sort.  split.  strings.tac.    tail.    tee
tr.    uniq.  wc.   
```

## 8. Traps


## 9. Functions
- sometimes shell scripts get complicated
- to avoid large and complicated scripts use functions
- you divide large scripts into small chunks/entities called **functions**
- functions make shell script modular and easy to use
- functions avoid repetitive code
- function performs a specific task
- function used like normal command
- in other high level programming languages function is also known as precedure, method, subrotine or routine

One line functions inside { ... } must end with a semicolon. Otherwise you get an error on screen:
```
xrpm() { rpm2cpio "$1" | cpio -idvm; }
```
#### writing functions
```
name() {
  command list;
}
```
- you can create a function file
- all shell functions are tested as a command
- you must define a function at the start of a script
- you must load a function file at the start of a script using *source* command
- you can all funciton like normal command

#### pass arguments into a function
- shell functions have their own command line argument
- use variable $1, $2, ..., $n to access argument passed to the function

```
name() {
  arg1=$1
  arg2=$2
  comand on $arg1
}
```

#### local variable
```
local var=value
```
#### Shell functions library
- you can store all your function in a function files called function library
- you can load all function into the current script or te command promt
- the syntax is as follows to load all functions

```
. /path/to/your/functions.sh
```
example
```
#!/bin/bash
# set variables
declare -r TRUE=0
declare -r FALSE=1
declare -r PASSWD_FILE=/etc/passwd
#######################
function to_lower()
{
      local str="$@"
      local output
      output=$(tr '[A-Z]' '[a-z]'<<<"${str}")
      echo $output
}
```

how do i load myfunction.sh into the script
```
#!/bin/bash
# Load the myfunctions.sh
# My local path is /home/vivek/lsst2/myfunctions.sh
. /home/vivek/lsst2/myfunctions.sh
```
#### source command
- the source command can be used to load any functions file into the current shell script or command prompt
- it read and execute commands from given FILENAME and return
- the pathnames in $PATH are used to find the directory containing FILENAME

## 10. Interactive Scripts
#### Menu driven scripts

#### bash display dialog boxes
```
brew install dialog
```
## 11. Shell scripting help


## 12. Other information

```
#!/bin/bash
set -e #设置该选项后，当脚本中任何以一个命令执行返回的状态码不为0时就退出整个脚本
set -u #设置该选项后，当脚本在执行过程中尝试使用未定义过的变量时，报错并退出运行整个脚本
```
单分支if判断
```
if [ condition ]
then
  command
fi
```
双分枝if判断
```
if [ condition ]
then
  command
else
  command
fi
```
多分枝if判断
```
if [ condition ]
then
  command
elif [ condition ]
then
  command
...
fi
```

#### processing files with bash using for loops and globbing
There are three essential parts to creating a pipeline to process a set of files:
1. selecting which files to apply the commands to
2. looping over the data and applying the commands
3. keeping track of the names of anyout put files created

There are different computational tricks to achieve each of these tasks. Let's first look at the simple ways to select which files to apply commands to:
- approaches that start with file comtaining information about samples
- approaches that select files in directories using some criteria



