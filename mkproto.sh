protoc -I=./schemes --cpp_out=./src/schemes --python_out=./src/schemes ./schemes/ccset.proto
protoc -I=./schemes --cpp_out=./src/schemes --python_out=./src/schemes ./schemes/histogramset.proto
capnp compile -oc++:src schemes/tags.capnp
