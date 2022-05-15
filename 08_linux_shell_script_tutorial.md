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
- 

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


## 4. Conditionals exections(Decision making)


## 5. Bash Loops


## 6. Shell redirection


## 7. Pipes and filters


## 8. Traps


## 9. Functions


## 10. Interactive Scripts


## 11. Shell scripting help
