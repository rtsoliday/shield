#!/bin/sh  
# \
exec tclsh "$0" "$@"

set version 1.1
set name shield-$version
puts "Building $name RPM"

exec ./rpmdev-setuptree
exec cp -f shield.spec $env(HOME)/rpmbuild/SPECS/
exec rm -rf $env(HOME)/rpmbuild/BUILD/$name
exec mkdir $env(HOME)/rpmbuild/BUILD/$name
set binFiles "[glob ../bin/Linux-x86_64/shield]"
foreach f $binFiles {
  exec chmod a+rx $f
  exec chmod a-w $f
  exec cp -f $f $env(HOME)/rpmbuild/BUILD/${name}/
}
exec cp -f ../src/shieldData.sdds $env(HOME)/rpmbuild/BUILD/${name}/
cd $env(HOME)/rpmbuild/BUILD
exec tar -cvf ../SOURCES/${name}.tar $name
exec rm -f ../SOURCES/${name}.tar.gz
exec gzip -9 ../SOURCES/${name}.tar
cd ../SPECS
if {[catch {exec rpmbuild -bb --quiet --clean --target x86_64 \
                 --buildroot $env(HOME)/rpmbuild/BUILDROOT shield.spec} results]} {
}
exec rm -f ../SOURCES/${name}.tar.gz
puts $results


puts "New RPM file in ~/rpmbuild/RPMS/x86_64"






