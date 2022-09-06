# -*- coding: utf-8 -*-
# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: schemes/ccset.proto

from google.protobuf import descriptor as _descriptor
from google.protobuf import message as _message
from google.protobuf import reflection as _reflection
from google.protobuf import symbol_database as _symbol_database
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()




DESCRIPTOR = _descriptor.FileDescriptor(
  name='schemes/ccset.proto',
  package='ccset',
  syntax='proto3',
  serialized_options=None,
  create_key=_descriptor._internal_create_key,
  serialized_pb=b'\n\x13schemes/ccset.proto\x12\x05\x63\x63set\"\xe7\x01\n\nccset_data\x12\x0f\n\x07offsets\x18\x01 \x03(\x01\x12\x0c\n\x04\x63\x63_h\x18\x02 \x03(\x03\x12\x0c\n\x04\x63\x63_v\x18\x03 \x03(\x03\x12)\n\tcc_h_tags\x18\x04 \x03(\x0b\x32\x16.ccset.ccset_data.darr\x12)\n\tcc_v_tags\x18\x05 \x03(\x0b\x32\x16.ccset.ccset_data.darr\x12\x10\n\x08meastime\x18\x06 \x01(\x01\x12\x18\n\x10num_trigger_tags\x18\x07 \x01(\x03\x12\x15\n\rnum_fpga_tags\x18\x08 \x01(\x03\x1a\x13\n\x04\x64\x61rr\x12\x0b\n\x03\x61rr\x18\x01 \x03(\x01\x62\x06proto3'
)




_CCSET_DATA_DARR = _descriptor.Descriptor(
  name='darr',
  full_name='ccset.ccset_data.darr',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  create_key=_descriptor._internal_create_key,
  fields=[
    _descriptor.FieldDescriptor(
      name='arr', full_name='ccset.ccset_data.darr.arr', index=0,
      number=1, type=1, cpp_type=5, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
  ],
  extensions=[
  ],
  nested_types=[],
  enum_types=[
  ],
  serialized_options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=243,
  serialized_end=262,
)

_CCSET_DATA = _descriptor.Descriptor(
  name='ccset_data',
  full_name='ccset.ccset_data',
  filename=None,
  file=DESCRIPTOR,
  containing_type=None,
  create_key=_descriptor._internal_create_key,
  fields=[
    _descriptor.FieldDescriptor(
      name='offsets', full_name='ccset.ccset_data.offsets', index=0,
      number=1, type=1, cpp_type=5, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='cc_h', full_name='ccset.ccset_data.cc_h', index=1,
      number=2, type=3, cpp_type=2, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='cc_v', full_name='ccset.ccset_data.cc_v', index=2,
      number=3, type=3, cpp_type=2, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='cc_h_tags', full_name='ccset.ccset_data.cc_h_tags', index=3,
      number=4, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='cc_v_tags', full_name='ccset.ccset_data.cc_v_tags', index=4,
      number=5, type=11, cpp_type=10, label=3,
      has_default_value=False, default_value=[],
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='meastime', full_name='ccset.ccset_data.meastime', index=5,
      number=6, type=1, cpp_type=5, label=1,
      has_default_value=False, default_value=float(0),
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='num_trigger_tags', full_name='ccset.ccset_data.num_trigger_tags', index=6,
      number=7, type=3, cpp_type=2, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
    _descriptor.FieldDescriptor(
      name='num_fpga_tags', full_name='ccset.ccset_data.num_fpga_tags', index=7,
      number=8, type=3, cpp_type=2, label=1,
      has_default_value=False, default_value=0,
      message_type=None, enum_type=None, containing_type=None,
      is_extension=False, extension_scope=None,
      serialized_options=None, file=DESCRIPTOR,  create_key=_descriptor._internal_create_key),
  ],
  extensions=[
  ],
  nested_types=[_CCSET_DATA_DARR, ],
  enum_types=[
  ],
  serialized_options=None,
  is_extendable=False,
  syntax='proto3',
  extension_ranges=[],
  oneofs=[
  ],
  serialized_start=31,
  serialized_end=262,
)

_CCSET_DATA_DARR.containing_type = _CCSET_DATA
_CCSET_DATA.fields_by_name['cc_h_tags'].message_type = _CCSET_DATA_DARR
_CCSET_DATA.fields_by_name['cc_v_tags'].message_type = _CCSET_DATA_DARR
DESCRIPTOR.message_types_by_name['ccset_data'] = _CCSET_DATA
_sym_db.RegisterFileDescriptor(DESCRIPTOR)

ccset_data = _reflection.GeneratedProtocolMessageType('ccset_data', (_message.Message,), {

  'darr' : _reflection.GeneratedProtocolMessageType('darr', (_message.Message,), {
    'DESCRIPTOR' : _CCSET_DATA_DARR,
    '__module__' : 'schemes.ccset_pb2'
    # @@protoc_insertion_point(class_scope:ccset.ccset_data.darr)
    })
  ,
  'DESCRIPTOR' : _CCSET_DATA,
  '__module__' : 'schemes.ccset_pb2'
  # @@protoc_insertion_point(class_scope:ccset.ccset_data)
  })
_sym_db.RegisterMessage(ccset_data)
_sym_db.RegisterMessage(ccset_data.darr)


# @@protoc_insertion_point(module_scope)