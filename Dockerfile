FROM erdc/proteus:latest

MAINTAINER Proteus Developers <proteus@groups.google.com>

USER $NB_USER

COPY . $HOME
