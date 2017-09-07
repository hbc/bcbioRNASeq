#!/bin/sh

setup_git() {
  git config --global user.email "travis@travis-ci.org"
  git config --global user.name "Travis CI"
}

commit_website_files() {
  git fetch origin gh-pages
  git checkout gh-pages
  git pull
  cp report/*.html .
  git add *.html
  git commit --message "Travis build: $TRAVIS_BUILD_NUMBER"
}

upload_files() {
  git remote add origin-pages https://${GITHUB_TOKE}@github.com/bcbio/bcbio_rnaseq_output_example.git > /dev/null 2>&1
  git push --force --quiet --set-upstream origin-pages gh-pages 
}

setup_git
commit_website_files
upload_files
