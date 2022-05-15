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


## 6. Shell redirection


## 7. Pipes and filters


## 8. Traps


## 9. Functions


## 10. Interactive Scripts


## 11. Shell scripting help