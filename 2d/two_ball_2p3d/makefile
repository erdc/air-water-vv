a:
	python setup.py build_ext -i
	otool -L Chrono3D.so
b:
	install_name_tool -change libChronoEngine.dylib @rpath/libChronoEngine.dylib Chrono3D.so
	otool -L Chrono3D.so
