FROM debian:9

# Set the working directory.
WORKDIR /usr/src/app

### dependencies should almost be all installed
COPY cpp_used .

### install dependencies
RUN apt-get update && \
    apt-get -y install gcc-6 && \
    apt-get -y install g++-6
RUN apt-get -y install gawk

### compile asilo
RUN gcc-6 -O3 -lstdc++ --std=c++11 -g -lm \
        asilo.1.0.cpp\
        dependencies/readGenelist.cpp \
        dependencies/filterFunctions.cpp \
        dependencies/atavCategories.cpp \
        dependencies/basicFunctions.cpp \
        dependencies/arrangeResults.cpp \
        dependencies/trioComphet.cpp \
        -o asilo.1.0

### compile mergeGenelists
RUN gcc-6 -O3 -lstdc++ --std=c++11 -g \
        mergeGenelists.0.cpp\
        -o mergeGenelists.0

### copy AWK script for VCF2CSV conversion
COPY header-converter.awk .

### copy and start waiter script
COPY waiter-launcher-v0.sh .

### select ENTRYPOINT
ENTRYPOINT bash /usr/src/app/waiter-launcher-v0.sh
