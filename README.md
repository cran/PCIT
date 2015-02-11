[![Build Status](https://travis-ci.org/nathanhaigh/pcit.png?branch=master)](https://travis-ci.org/nathanhaigh/pcit)
pcit
====

PCIT: an R package for weighted gene co-expression networks based on partial correlation and information theory approaches

Building and Submitting Packages
================================
 1) Update "Version" and "Date" in the DESCRIPTION file
 2) Add info about changes since last release to ChangeLog file
 3) `R CMD build pcit`
 4) Check the built tarball: `R CMD check --as-cran PCIT_*.tar.gz`
 5) Submit the tarball at http://cran.r-project.org/submit.html or via FTP uisng instructions here: http://cran.r-project.org/banner.shtml#submitting
 6) Tag the release on github:

```
# Parse out the version info from the DESCRIOTION file
version=$(grep "^Version:" DESCRIPTION | awk '{print $2}')
# Tag the current branch using the Changelog entry for this version, asthe message
git tag -a "v${version}" -F <(perl -e '$v_encountered=0; while(<>){if (/'${version}'/){$v_encountered=1;print;next}; last if /^\d/; print}' ChangeLog)
# Push the tag up to github
git push origin "v${version}"

```

