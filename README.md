# vk_pygeom
This is a collection of geometry related python packages

# License
This software is released under the condition described in the LICENSE.

## install
```
python3 -m pip install vk2gpz-pygeom
```
or
```
pip install git+ssh://git@github.com/vk2gpz/vk_pygeom.git
```
(notice that ':' in "pip install git+ssh://git@github.com:vk2gpz/vk_pygeom.git" needs to be changed to '/')


## test
python
```
>>> import vk2gpz.geom.grid.geodesicdome as gd
>>> gd.main()
```

## build distribution archive
```
python3 -m pip install --upgrade build
python3 -m build
```

## uploading the distribution archive (change testpypi to pypi for official upload)
```
python3 -m pip install --upgrade twine
python3 -m twine upload --repository testpypi dist/*
```

## install pacakge from the project locally
```angular2html
pip install dynamodb_dataframes --no-index --find-links ./dist
```