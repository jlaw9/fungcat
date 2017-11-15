# utils
A collection of utility functions used by several other repositories and software packages.

## git subtree
To add this repository as a `git subtree` in the `src/utils` folder of another git project repository:

``` sh
git subtree add --prefix src/utils https://github.com/Murali-group/utils.git master --squash
```

To pull changes to the subtree:

``` sh
git subtree pull --prefix src/utils https://github.com/Murali-group/utils.git master --squash
```

## wiki
I started a wiki to describe the utils we add here. Feel free to add to it.

- https://github.com/Murali-group/utils/wiki/Baobab-info
