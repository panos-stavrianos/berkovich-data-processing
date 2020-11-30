# Installation for windows

## Install Python
Download from the original [site](https://www.python.org/downloads/release/python-390)
 and install python

## Install PIP
Either download pip from this [site](https://bootstrap.pypa.io/get-pip.py) (right click, save as)

 or use curl in command line
```
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
```

then run the file with
```
python get-pip.py
```
or with
```
py get-pip.py
```
## Download the project
To download the project either use the `Download zip` from github or

open terminal and navigate to the path you want and use git to download it 
(you must have [git](https://git-scm.com/download/win) installed)

```
cd C:\Users\some_user\apps
git clone https://github.com/panos-stavrianos/berkovich-data-processing
```

## Installing virtualenv
Open terminal and navigate to the project path
```
cd C:\Users\some_user\apps\berkovich-data-processing
python3 -m venv venv
```

## Install requirements.txt
```
.\venv\Scripts\activate
pip install -r requirements.txt
deactivate
```

## Start the app




