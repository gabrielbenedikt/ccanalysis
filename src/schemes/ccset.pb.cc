// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: schemes/ccset.proto

#include "schemes/ccset.pb.h"

#include <algorithm>

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/extension_set.h>
#include <google/protobuf/wire_format_lite.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/generated_message_reflection.h>
#include <google/protobuf/reflection_ops.h>
#include <google/protobuf/wire_format.h>
// @@protoc_insertion_point(includes)
#include <google/protobuf/port_def.inc>
extern PROTOBUF_INTERNAL_EXPORT_schemes_2fccset_2eproto ::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_ccset_data_darr_schemes_2fccset_2eproto;
namespace ccset {
class ccset_data_darrDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<ccset_data_darr> _instance;
} _ccset_data_darr_default_instance_;
class ccset_dataDefaultTypeInternal {
 public:
  ::PROTOBUF_NAMESPACE_ID::internal::ExplicitlyConstructed<ccset_data> _instance;
} _ccset_data_default_instance_;
}  // namespace ccset
static void InitDefaultsscc_info_ccset_data_schemes_2fccset_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::ccset::_ccset_data_default_instance_;
    new (ptr) ::ccset::ccset_data();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::ccset::ccset_data::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<1> scc_info_ccset_data_schemes_2fccset_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 1, 0, InitDefaultsscc_info_ccset_data_schemes_2fccset_2eproto}, {
      &scc_info_ccset_data_darr_schemes_2fccset_2eproto.base,}};

static void InitDefaultsscc_info_ccset_data_darr_schemes_2fccset_2eproto() {
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  {
    void* ptr = &::ccset::_ccset_data_darr_default_instance_;
    new (ptr) ::ccset::ccset_data_darr();
    ::PROTOBUF_NAMESPACE_ID::internal::OnShutdownDestroyMessage(ptr);
  }
  ::ccset::ccset_data_darr::InitAsDefaultInstance();
}

::PROTOBUF_NAMESPACE_ID::internal::SCCInfo<0> scc_info_ccset_data_darr_schemes_2fccset_2eproto =
    {{ATOMIC_VAR_INIT(::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase::kUninitialized), 0, 0, InitDefaultsscc_info_ccset_data_darr_schemes_2fccset_2eproto}, {}};

static ::PROTOBUF_NAMESPACE_ID::Metadata file_level_metadata_schemes_2fccset_2eproto[2];
static constexpr ::PROTOBUF_NAMESPACE_ID::EnumDescriptor const** file_level_enum_descriptors_schemes_2fccset_2eproto = nullptr;
static constexpr ::PROTOBUF_NAMESPACE_ID::ServiceDescriptor const** file_level_service_descriptors_schemes_2fccset_2eproto = nullptr;

const ::PROTOBUF_NAMESPACE_ID::uint32 TableStruct_schemes_2fccset_2eproto::offsets[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  ~0u,  // no _has_bits_
  PROTOBUF_FIELD_OFFSET(::ccset::ccset_data_darr, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::ccset::ccset_data_darr, arr_),
  ~0u,  // no _has_bits_
  PROTOBUF_FIELD_OFFSET(::ccset::ccset_data, _internal_metadata_),
  ~0u,  // no _extensions_
  ~0u,  // no _oneof_case_
  ~0u,  // no _weak_field_map_
  PROTOBUF_FIELD_OFFSET(::ccset::ccset_data, offsets_),
  PROTOBUF_FIELD_OFFSET(::ccset::ccset_data, cc_h_),
  PROTOBUF_FIELD_OFFSET(::ccset::ccset_data, cc_v_),
  PROTOBUF_FIELD_OFFSET(::ccset::ccset_data, cc_h_tags_),
  PROTOBUF_FIELD_OFFSET(::ccset::ccset_data, cc_v_tags_),
  PROTOBUF_FIELD_OFFSET(::ccset::ccset_data, meastime_),
  PROTOBUF_FIELD_OFFSET(::ccset::ccset_data, num_trigger_tags_),
  PROTOBUF_FIELD_OFFSET(::ccset::ccset_data, num_fpga_tags_),
};
static const ::PROTOBUF_NAMESPACE_ID::internal::MigrationSchema schemas[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) = {
  { 0, -1, sizeof(::ccset::ccset_data_darr)},
  { 6, -1, sizeof(::ccset::ccset_data)},
};

static ::PROTOBUF_NAMESPACE_ID::Message const * const file_default_instances[] = {
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::ccset::_ccset_data_darr_default_instance_),
  reinterpret_cast<const ::PROTOBUF_NAMESPACE_ID::Message*>(&::ccset::_ccset_data_default_instance_),
};

const char descriptor_table_protodef_schemes_2fccset_2eproto[] PROTOBUF_SECTION_VARIABLE(protodesc_cold) =
  "\n\023schemes/ccset.proto\022\005ccset\"\347\001\n\nccset_d"
  "ata\022\017\n\007offsets\030\001 \003(\001\022\014\n\004cc_h\030\002 \003(\003\022\014\n\004cc"
  "_v\030\003 \003(\003\022)\n\tcc_h_tags\030\004 \003(\0132\026.ccset.ccse"
  "t_data.darr\022)\n\tcc_v_tags\030\005 \003(\0132\026.ccset.c"
  "cset_data.darr\022\020\n\010meastime\030\006 \001(\001\022\030\n\020num_"
  "trigger_tags\030\007 \001(\003\022\025\n\rnum_fpga_tags\030\010 \001("
  "\003\032\023\n\004darr\022\013\n\003arr\030\001 \003(\001b\006proto3"
  ;
static const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable*const descriptor_table_schemes_2fccset_2eproto_deps[1] = {
};
static ::PROTOBUF_NAMESPACE_ID::internal::SCCInfoBase*const descriptor_table_schemes_2fccset_2eproto_sccs[2] = {
  &scc_info_ccset_data_schemes_2fccset_2eproto.base,
  &scc_info_ccset_data_darr_schemes_2fccset_2eproto.base,
};
static ::PROTOBUF_NAMESPACE_ID::internal::once_flag descriptor_table_schemes_2fccset_2eproto_once;
const ::PROTOBUF_NAMESPACE_ID::internal::DescriptorTable descriptor_table_schemes_2fccset_2eproto = {
  false, false, descriptor_table_protodef_schemes_2fccset_2eproto, "schemes/ccset.proto", 270,
  &descriptor_table_schemes_2fccset_2eproto_once, descriptor_table_schemes_2fccset_2eproto_sccs, descriptor_table_schemes_2fccset_2eproto_deps, 2, 0,
  schemas, file_default_instances, TableStruct_schemes_2fccset_2eproto::offsets,
  file_level_metadata_schemes_2fccset_2eproto, 2, file_level_enum_descriptors_schemes_2fccset_2eproto, file_level_service_descriptors_schemes_2fccset_2eproto,
};

// Force running AddDescriptors() at dynamic initialization time.
static bool dynamic_init_dummy_schemes_2fccset_2eproto = (static_cast<void>(::PROTOBUF_NAMESPACE_ID::internal::AddDescriptors(&descriptor_table_schemes_2fccset_2eproto)), true);
namespace ccset {

// ===================================================================

void ccset_data_darr::InitAsDefaultInstance() {
}
class ccset_data_darr::_Internal {
 public:
};

ccset_data_darr::ccset_data_darr(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  arr_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:ccset.ccset_data.darr)
}
ccset_data_darr::ccset_data_darr(const ccset_data_darr& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      arr_(from.arr_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  // @@protoc_insertion_point(copy_constructor:ccset.ccset_data.darr)
}

void ccset_data_darr::SharedCtor() {
}

ccset_data_darr::~ccset_data_darr() {
  // @@protoc_insertion_point(destructor:ccset.ccset_data.darr)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void ccset_data_darr::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void ccset_data_darr::ArenaDtor(void* object) {
  ccset_data_darr* _this = reinterpret_cast< ccset_data_darr* >(object);
  (void)_this;
}
void ccset_data_darr::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void ccset_data_darr::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const ccset_data_darr& ccset_data_darr::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_ccset_data_darr_schemes_2fccset_2eproto.base);
  return *internal_default_instance();
}


void ccset_data_darr::Clear() {
// @@protoc_insertion_point(message_clear_start:ccset.ccset_data.darr)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  arr_.Clear();
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* ccset_data_darr::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // repeated double arr = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 10)) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedDoubleParser(_internal_mutable_arr(), ptr, ctx);
          CHK_(ptr);
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 9) {
          _internal_add_arr(::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr));
          ptr += sizeof(double);
        } else goto handle_unusual;
        continue;
      default: {
      handle_unusual:
        if ((tag & 7) == 4 || tag == 0) {
          ctx->SetLastTag(tag);
          goto success;
        }
        ptr = UnknownFieldParse(tag,
            _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(),
            ptr, ctx);
        CHK_(ptr != nullptr);
        continue;
      }
    }  // switch
  }  // while
success:
  return ptr;
failure:
  ptr = nullptr;
  goto success;
#undef CHK_
}

::PROTOBUF_NAMESPACE_ID::uint8* ccset_data_darr::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:ccset.ccset_data.darr)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // repeated double arr = 1;
  if (this->_internal_arr_size() > 0) {
    target = stream->WriteFixedPacked(1, _internal_arr(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:ccset.ccset_data.darr)
  return target;
}

size_t ccset_data_darr::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:ccset.ccset_data.darr)
  size_t total_size = 0;

  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated double arr = 1;
  {
    unsigned int count = static_cast<unsigned int>(this->_internal_arr_size());
    size_t data_size = 8UL * count;
    if (data_size > 0) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::Int32Size(
            static_cast<::PROTOBUF_NAMESPACE_ID::int32>(data_size));
    }
    int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(data_size);
    _arr_cached_byte_size_.store(cached_size,
                                    std::memory_order_relaxed);
    total_size += data_size;
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    return ::PROTOBUF_NAMESPACE_ID::internal::ComputeUnknownFieldsSize(
        _internal_metadata_, total_size, &_cached_size_);
  }
  int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void ccset_data_darr::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:ccset.ccset_data.darr)
  GOOGLE_DCHECK_NE(&from, this);
  const ccset_data_darr* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<ccset_data_darr>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:ccset.ccset_data.darr)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:ccset.ccset_data.darr)
    MergeFrom(*source);
  }
}

void ccset_data_darr::MergeFrom(const ccset_data_darr& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:ccset.ccset_data.darr)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  arr_.MergeFrom(from.arr_);
}

void ccset_data_darr::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:ccset.ccset_data.darr)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void ccset_data_darr::CopyFrom(const ccset_data_darr& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:ccset.ccset_data.darr)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool ccset_data_darr::IsInitialized() const {
  return true;
}

void ccset_data_darr::InternalSwap(ccset_data_darr* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  arr_.InternalSwap(&other->arr_);
}

::PROTOBUF_NAMESPACE_ID::Metadata ccset_data_darr::GetMetadata() const {
  return GetMetadataStatic();
}


// ===================================================================

void ccset_data::InitAsDefaultInstance() {
}
class ccset_data::_Internal {
 public:
};

ccset_data::ccset_data(::PROTOBUF_NAMESPACE_ID::Arena* arena)
  : ::PROTOBUF_NAMESPACE_ID::Message(arena),
  offsets_(arena),
  cc_h_(arena),
  cc_v_(arena),
  cc_h_tags_(arena),
  cc_v_tags_(arena) {
  SharedCtor();
  RegisterArenaDtor(arena);
  // @@protoc_insertion_point(arena_constructor:ccset.ccset_data)
}
ccset_data::ccset_data(const ccset_data& from)
  : ::PROTOBUF_NAMESPACE_ID::Message(),
      offsets_(from.offsets_),
      cc_h_(from.cc_h_),
      cc_v_(from.cc_v_),
      cc_h_tags_(from.cc_h_tags_),
      cc_v_tags_(from.cc_v_tags_) {
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::memcpy(&meastime_, &from.meastime_,
    static_cast<size_t>(reinterpret_cast<char*>(&num_fpga_tags_) -
    reinterpret_cast<char*>(&meastime_)) + sizeof(num_fpga_tags_));
  // @@protoc_insertion_point(copy_constructor:ccset.ccset_data)
}

void ccset_data::SharedCtor() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&scc_info_ccset_data_schemes_2fccset_2eproto.base);
  ::memset(&meastime_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&num_fpga_tags_) -
      reinterpret_cast<char*>(&meastime_)) + sizeof(num_fpga_tags_));
}

ccset_data::~ccset_data() {
  // @@protoc_insertion_point(destructor:ccset.ccset_data)
  SharedDtor();
  _internal_metadata_.Delete<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

void ccset_data::SharedDtor() {
  GOOGLE_DCHECK(GetArena() == nullptr);
}

void ccset_data::ArenaDtor(void* object) {
  ccset_data* _this = reinterpret_cast< ccset_data* >(object);
  (void)_this;
}
void ccset_data::RegisterArenaDtor(::PROTOBUF_NAMESPACE_ID::Arena*) {
}
void ccset_data::SetCachedSize(int size) const {
  _cached_size_.Set(size);
}
const ccset_data& ccset_data::default_instance() {
  ::PROTOBUF_NAMESPACE_ID::internal::InitSCC(&::scc_info_ccset_data_schemes_2fccset_2eproto.base);
  return *internal_default_instance();
}


void ccset_data::Clear() {
// @@protoc_insertion_point(message_clear_start:ccset.ccset_data)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  offsets_.Clear();
  cc_h_.Clear();
  cc_v_.Clear();
  cc_h_tags_.Clear();
  cc_v_tags_.Clear();
  ::memset(&meastime_, 0, static_cast<size_t>(
      reinterpret_cast<char*>(&num_fpga_tags_) -
      reinterpret_cast<char*>(&meastime_)) + sizeof(num_fpga_tags_));
  _internal_metadata_.Clear<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>();
}

const char* ccset_data::_InternalParse(const char* ptr, ::PROTOBUF_NAMESPACE_ID::internal::ParseContext* ctx) {
#define CHK_(x) if (PROTOBUF_PREDICT_FALSE(!(x))) goto failure
  ::PROTOBUF_NAMESPACE_ID::Arena* arena = GetArena(); (void)arena;
  while (!ctx->Done(&ptr)) {
    ::PROTOBUF_NAMESPACE_ID::uint32 tag;
    ptr = ::PROTOBUF_NAMESPACE_ID::internal::ReadTag(ptr, &tag);
    CHK_(ptr);
    switch (tag >> 3) {
      // repeated double offsets = 1;
      case 1:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 10)) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedDoubleParser(_internal_mutable_offsets(), ptr, ctx);
          CHK_(ptr);
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 9) {
          _internal_add_offsets(::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr));
          ptr += sizeof(double);
        } else goto handle_unusual;
        continue;
      // repeated int64 cc_h = 2;
      case 2:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 18)) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedInt64Parser(_internal_mutable_cc_h(), ptr, ctx);
          CHK_(ptr);
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 16) {
          _internal_add_cc_h(::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr));
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // repeated int64 cc_v = 3;
      case 3:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 26)) {
          ptr = ::PROTOBUF_NAMESPACE_ID::internal::PackedInt64Parser(_internal_mutable_cc_v(), ptr, ctx);
          CHK_(ptr);
        } else if (static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 24) {
          _internal_add_cc_v(::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr));
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // repeated .ccset.ccset_data.darr cc_h_tags = 4;
      case 4:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 34)) {
          ptr -= 1;
          do {
            ptr += 1;
            ptr = ctx->ParseMessage(_internal_add_cc_h_tags(), ptr);
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<34>(ptr));
        } else goto handle_unusual;
        continue;
      // repeated .ccset.ccset_data.darr cc_v_tags = 5;
      case 5:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 42)) {
          ptr -= 1;
          do {
            ptr += 1;
            ptr = ctx->ParseMessage(_internal_add_cc_v_tags(), ptr);
            CHK_(ptr);
            if (!ctx->DataAvailable(ptr)) break;
          } while (::PROTOBUF_NAMESPACE_ID::internal::ExpectTag<42>(ptr));
        } else goto handle_unusual;
        continue;
      // double meastime = 6;
      case 6:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 49)) {
          meastime_ = ::PROTOBUF_NAMESPACE_ID::internal::UnalignedLoad<double>(ptr);
          ptr += sizeof(double);
        } else goto handle_unusual;
        continue;
      // int64 num_trigger_tags = 7;
      case 7:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 56)) {
          num_trigger_tags_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      // int64 num_fpga_tags = 8;
      case 8:
        if (PROTOBUF_PREDICT_TRUE(static_cast<::PROTOBUF_NAMESPACE_ID::uint8>(tag) == 64)) {
          num_fpga_tags_ = ::PROTOBUF_NAMESPACE_ID::internal::ReadVarint64(&ptr);
          CHK_(ptr);
        } else goto handle_unusual;
        continue;
      default: {
      handle_unusual:
        if ((tag & 7) == 4 || tag == 0) {
          ctx->SetLastTag(tag);
          goto success;
        }
        ptr = UnknownFieldParse(tag,
            _internal_metadata_.mutable_unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(),
            ptr, ctx);
        CHK_(ptr != nullptr);
        continue;
      }
    }  // switch
  }  // while
success:
  return ptr;
failure:
  ptr = nullptr;
  goto success;
#undef CHK_
}

::PROTOBUF_NAMESPACE_ID::uint8* ccset_data::_InternalSerialize(
    ::PROTOBUF_NAMESPACE_ID::uint8* target, ::PROTOBUF_NAMESPACE_ID::io::EpsCopyOutputStream* stream) const {
  // @@protoc_insertion_point(serialize_to_array_start:ccset.ccset_data)
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  // repeated double offsets = 1;
  if (this->_internal_offsets_size() > 0) {
    target = stream->WriteFixedPacked(1, _internal_offsets(), target);
  }

  // repeated int64 cc_h = 2;
  {
    int byte_size = _cc_h_cached_byte_size_.load(std::memory_order_relaxed);
    if (byte_size > 0) {
      target = stream->WriteInt64Packed(
          2, _internal_cc_h(), byte_size, target);
    }
  }

  // repeated int64 cc_v = 3;
  {
    int byte_size = _cc_v_cached_byte_size_.load(std::memory_order_relaxed);
    if (byte_size > 0) {
      target = stream->WriteInt64Packed(
          3, _internal_cc_v(), byte_size, target);
    }
  }

  // repeated .ccset.ccset_data.darr cc_h_tags = 4;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_cc_h_tags_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(4, this->_internal_cc_h_tags(i), target, stream);
  }

  // repeated .ccset.ccset_data.darr cc_v_tags = 5;
  for (unsigned int i = 0,
      n = static_cast<unsigned int>(this->_internal_cc_v_tags_size()); i < n; i++) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      InternalWriteMessage(5, this->_internal_cc_v_tags(i), target, stream);
  }

  // double meastime = 6;
  if (!(this->meastime() <= 0 && this->meastime() >= 0)) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteDoubleToArray(6, this->_internal_meastime(), target);
  }

  // int64 num_trigger_tags = 7;
  if (this->num_trigger_tags() != 0) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteInt64ToArray(7, this->_internal_num_trigger_tags(), target);
  }

  // int64 num_fpga_tags = 8;
  if (this->num_fpga_tags() != 0) {
    target = stream->EnsureSpace(target);
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::WriteInt64ToArray(8, this->_internal_num_fpga_tags(), target);
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    target = ::PROTOBUF_NAMESPACE_ID::internal::WireFormat::InternalSerializeUnknownFieldsToArray(
        _internal_metadata_.unknown_fields<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(::PROTOBUF_NAMESPACE_ID::UnknownFieldSet::default_instance), target, stream);
  }
  // @@protoc_insertion_point(serialize_to_array_end:ccset.ccset_data)
  return target;
}

size_t ccset_data::ByteSizeLong() const {
// @@protoc_insertion_point(message_byte_size_start:ccset.ccset_data)
  size_t total_size = 0;

  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  // Prevent compiler warnings about cached_has_bits being unused
  (void) cached_has_bits;

  // repeated double offsets = 1;
  {
    unsigned int count = static_cast<unsigned int>(this->_internal_offsets_size());
    size_t data_size = 8UL * count;
    if (data_size > 0) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::Int32Size(
            static_cast<::PROTOBUF_NAMESPACE_ID::int32>(data_size));
    }
    int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(data_size);
    _offsets_cached_byte_size_.store(cached_size,
                                    std::memory_order_relaxed);
    total_size += data_size;
  }

  // repeated int64 cc_h = 2;
  {
    size_t data_size = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      Int64Size(this->cc_h_);
    if (data_size > 0) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::Int32Size(
            static_cast<::PROTOBUF_NAMESPACE_ID::int32>(data_size));
    }
    int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(data_size);
    _cc_h_cached_byte_size_.store(cached_size,
                                    std::memory_order_relaxed);
    total_size += data_size;
  }

  // repeated int64 cc_v = 3;
  {
    size_t data_size = ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::
      Int64Size(this->cc_v_);
    if (data_size > 0) {
      total_size += 1 +
        ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::Int32Size(
            static_cast<::PROTOBUF_NAMESPACE_ID::int32>(data_size));
    }
    int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(data_size);
    _cc_v_cached_byte_size_.store(cached_size,
                                    std::memory_order_relaxed);
    total_size += data_size;
  }

  // repeated .ccset.ccset_data.darr cc_h_tags = 4;
  total_size += 1UL * this->_internal_cc_h_tags_size();
  for (const auto& msg : this->cc_h_tags_) {
    total_size +=
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(msg);
  }

  // repeated .ccset.ccset_data.darr cc_v_tags = 5;
  total_size += 1UL * this->_internal_cc_v_tags_size();
  for (const auto& msg : this->cc_v_tags_) {
    total_size +=
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::MessageSize(msg);
  }

  // double meastime = 6;
  if (!(this->meastime() <= 0 && this->meastime() >= 0)) {
    total_size += 1 + 8;
  }

  // int64 num_trigger_tags = 7;
  if (this->num_trigger_tags() != 0) {
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::Int64Size(
        this->_internal_num_trigger_tags());
  }

  // int64 num_fpga_tags = 8;
  if (this->num_fpga_tags() != 0) {
    total_size += 1 +
      ::PROTOBUF_NAMESPACE_ID::internal::WireFormatLite::Int64Size(
        this->_internal_num_fpga_tags());
  }

  if (PROTOBUF_PREDICT_FALSE(_internal_metadata_.have_unknown_fields())) {
    return ::PROTOBUF_NAMESPACE_ID::internal::ComputeUnknownFieldsSize(
        _internal_metadata_, total_size, &_cached_size_);
  }
  int cached_size = ::PROTOBUF_NAMESPACE_ID::internal::ToCachedSize(total_size);
  SetCachedSize(cached_size);
  return total_size;
}

void ccset_data::MergeFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_merge_from_start:ccset.ccset_data)
  GOOGLE_DCHECK_NE(&from, this);
  const ccset_data* source =
      ::PROTOBUF_NAMESPACE_ID::DynamicCastToGenerated<ccset_data>(
          &from);
  if (source == nullptr) {
  // @@protoc_insertion_point(generalized_merge_from_cast_fail:ccset.ccset_data)
    ::PROTOBUF_NAMESPACE_ID::internal::ReflectionOps::Merge(from, this);
  } else {
  // @@protoc_insertion_point(generalized_merge_from_cast_success:ccset.ccset_data)
    MergeFrom(*source);
  }
}

void ccset_data::MergeFrom(const ccset_data& from) {
// @@protoc_insertion_point(class_specific_merge_from_start:ccset.ccset_data)
  GOOGLE_DCHECK_NE(&from, this);
  _internal_metadata_.MergeFrom<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(from._internal_metadata_);
  ::PROTOBUF_NAMESPACE_ID::uint32 cached_has_bits = 0;
  (void) cached_has_bits;

  offsets_.MergeFrom(from.offsets_);
  cc_h_.MergeFrom(from.cc_h_);
  cc_v_.MergeFrom(from.cc_v_);
  cc_h_tags_.MergeFrom(from.cc_h_tags_);
  cc_v_tags_.MergeFrom(from.cc_v_tags_);
  if (!(from.meastime() <= 0 && from.meastime() >= 0)) {
    _internal_set_meastime(from._internal_meastime());
  }
  if (from.num_trigger_tags() != 0) {
    _internal_set_num_trigger_tags(from._internal_num_trigger_tags());
  }
  if (from.num_fpga_tags() != 0) {
    _internal_set_num_fpga_tags(from._internal_num_fpga_tags());
  }
}

void ccset_data::CopyFrom(const ::PROTOBUF_NAMESPACE_ID::Message& from) {
// @@protoc_insertion_point(generalized_copy_from_start:ccset.ccset_data)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void ccset_data::CopyFrom(const ccset_data& from) {
// @@protoc_insertion_point(class_specific_copy_from_start:ccset.ccset_data)
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool ccset_data::IsInitialized() const {
  return true;
}

void ccset_data::InternalSwap(ccset_data* other) {
  using std::swap;
  _internal_metadata_.Swap<::PROTOBUF_NAMESPACE_ID::UnknownFieldSet>(&other->_internal_metadata_);
  offsets_.InternalSwap(&other->offsets_);
  cc_h_.InternalSwap(&other->cc_h_);
  cc_v_.InternalSwap(&other->cc_v_);
  cc_h_tags_.InternalSwap(&other->cc_h_tags_);
  cc_v_tags_.InternalSwap(&other->cc_v_tags_);
  ::PROTOBUF_NAMESPACE_ID::internal::memswap<
      PROTOBUF_FIELD_OFFSET(ccset_data, num_fpga_tags_)
      + sizeof(ccset_data::num_fpga_tags_)
      - PROTOBUF_FIELD_OFFSET(ccset_data, meastime_)>(
          reinterpret_cast<char*>(&meastime_),
          reinterpret_cast<char*>(&other->meastime_));
}

::PROTOBUF_NAMESPACE_ID::Metadata ccset_data::GetMetadata() const {
  return GetMetadataStatic();
}


// @@protoc_insertion_point(namespace_scope)
}  // namespace ccset
PROTOBUF_NAMESPACE_OPEN
template<> PROTOBUF_NOINLINE ::ccset::ccset_data_darr* Arena::CreateMaybeMessage< ::ccset::ccset_data_darr >(Arena* arena) {
  return Arena::CreateMessageInternal< ::ccset::ccset_data_darr >(arena);
}
template<> PROTOBUF_NOINLINE ::ccset::ccset_data* Arena::CreateMaybeMessage< ::ccset::ccset_data >(Arena* arena) {
  return Arena::CreateMessageInternal< ::ccset::ccset_data >(arena);
}
PROTOBUF_NAMESPACE_CLOSE

// @@protoc_insertion_point(global_scope)
#include <google/protobuf/port_undef.inc>
