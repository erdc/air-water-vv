#!/usr/bin/env python
# -*- coding: utf-8 -*- #
from __future__ import unicode_literals

AUTHOR = u'ERDC Computational Mechanics'
SITENAME = u'Verification and Validation for Air/Water Flow Models'
SITEURL = ''

PATH = 'content'

TIMEZONE = 'Europe/Paris'

DEFAULT_LANG = u'en'

# Feed generation is usually not desired when developing
FEED_ALL_ATOM = None
CATEGORY_FEED_ATOM = None
TRANSLATION_FEED_ATOM = None
AUTHOR_FEED_ATOM = None
AUTHOR_FEED_RSS = None

# Blogroll
LINKS = (('Pelican', 'http://getpelican.com/'),
         ('Python.org', 'http://python.org/'),
         ('Jinja2', 'http://jinja.pocoo.org/'))

DEFAULT_PAGINATION = 10

STATIC_PATHS = ['images',]

THEME = '../stack.Darwin/share/pelican-bootstrap3'
BOOTSTRAP_THEME = 'lumen'

# Uncomment following line if you want document-relative URLs when developing
#RELATIVE_URLS = True
