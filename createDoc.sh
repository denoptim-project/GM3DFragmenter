#!/bin/bash
find ./src -type f -name "*.java" | xargs javadoc -d doc -classpath lib/cdk-1.4.17.jar
