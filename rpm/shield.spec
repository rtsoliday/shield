Summary:	Accelerator code
Name:		shield
License:	EPICS Open license http://www.aps.anl.gov/epics/license/open.php
Group:		Applications/Databases
URL:		http://www.aps.anl.gov/asd/oag/oaghome.shtml
Packager:	Robert Soliday <soliday@aps.anl.gov>
Prefix:		%{_bindir}
Autoreq:	0
Version:	1.1
Release:	1
Source:		shield-1.1.tar.gz


%define debug_package %{nil}
%undefine __check_files
%description
Binary package for Shield. Shield does shielding analyses around
a high energy accelerator. This version is based off of SHIELD11 
from SLAC.

%prep
%setup

%build
%install
mkdir -p %{buildroot}%{_bindir}
mkdir -p %{buildroot}%{_prefix}/local/oag/apps/configData/shield
install -s -m 755 shield %{buildroot}%{_bindir}/shield
install -m 666 shieldData.sdds %{buildroot}%{_prefix}/local/oag/apps/configData/shield/shieldData.sdds

%files

%{_bindir}/shield
%{_prefix}/local/oag/apps/configData/shield/shieldData.sdds
