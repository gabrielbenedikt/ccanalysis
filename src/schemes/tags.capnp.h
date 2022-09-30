// Generated by Cap'n Proto compiler, DO NOT EDIT
// source: tags.capnp

#pragma once

#include <capnp/generated-header-support.h>
#include <kj/windows-sanity.h>

#if CAPNP_VERSION != 9001
#error "Version mismatch between generated code and library headers.  You must use the same version of the Cap'n Proto compiler and library."
#endif


CAPNP_BEGIN_HEADER

namespace capnp {
namespace schemas {

CAPNP_DECLARE_SCHEMA(b1642a9902d01394);
CAPNP_DECLARE_SCHEMA(8995b3a3aece585b);

}  // namespace schemas
}  // namespace capnp


struct Tags {
  Tags() = delete;

  class Reader;
  class Builder;
  class Pipeline;
  struct Tag;

  struct _capnpPrivate {
    CAPNP_DECLARE_STRUCT_HEADER(b1642a9902d01394, 0, 1)
    #if !CAPNP_LITE
    static constexpr ::capnp::_::RawBrandedSchema const* brand() { return &schema->defaultBrand; }
    #endif  // !CAPNP_LITE
  };
};

struct Tags::Tag {
  Tag() = delete;

  class Reader;
  class Builder;
  class Pipeline;

  struct _capnpPrivate {
    CAPNP_DECLARE_STRUCT_HEADER(8995b3a3aece585b, 2, 0)
    #if !CAPNP_LITE
    static constexpr ::capnp::_::RawBrandedSchema const* brand() { return &schema->defaultBrand; }
    #endif  // !CAPNP_LITE
  };
};

// =======================================================================================

class Tags::Reader {
public:
  typedef Tags Reads;

  Reader() = default;
  inline explicit Reader(::capnp::_::StructReader base): _reader(base) {}

  inline ::capnp::MessageSize totalSize() const {
    return _reader.totalSize().asPublic();
  }

#if !CAPNP_LITE
  inline ::kj::StringTree toString() const {
    return ::capnp::_::structString(_reader, *_capnpPrivate::brand());
  }
#endif  // !CAPNP_LITE

  inline bool hasTags() const;
  inline  ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>::Reader getTags() const;

private:
  ::capnp::_::StructReader _reader;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::ToDynamic_;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::_::PointerHelpers;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::List;
  friend class ::capnp::MessageBuilder;
  friend class ::capnp::Orphanage;
};

class Tags::Builder {
public:
  typedef Tags Builds;

  Builder() = delete;  // Deleted to discourage incorrect usage.
                       // You can explicitly initialize to nullptr instead.
  inline Builder(decltype(nullptr)) {}
  inline explicit Builder(::capnp::_::StructBuilder base): _builder(base) {}
  inline operator Reader() const { return Reader(_builder.asReader()); }
  inline Reader asReader() const { return *this; }

  inline ::capnp::MessageSize totalSize() const { return asReader().totalSize(); }
#if !CAPNP_LITE
  inline ::kj::StringTree toString() const { return asReader().toString(); }
#endif  // !CAPNP_LITE

  inline bool hasTags();
  inline  ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>::Builder getTags();
  inline void setTags( ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>::Reader value);
  inline void setTags(::kj::ArrayPtr<const  ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>::Reader> value);
  inline  ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>::Builder initTags(unsigned int size);
  inline void adoptTags(::capnp::Orphan< ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>>&& value);
  inline ::capnp::Orphan< ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>> disownTags();

private:
  ::capnp::_::StructBuilder _builder;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::ToDynamic_;
  friend class ::capnp::Orphanage;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::_::PointerHelpers;
};

#if !CAPNP_LITE
class Tags::Pipeline {
public:
  typedef Tags Pipelines;

  inline Pipeline(decltype(nullptr)): _typeless(nullptr) {}
  inline explicit Pipeline(::capnp::AnyPointer::Pipeline&& typeless)
      : _typeless(kj::mv(typeless)) {}

private:
  ::capnp::AnyPointer::Pipeline _typeless;
  friend class ::capnp::PipelineHook;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::ToDynamic_;
};
#endif  // !CAPNP_LITE

class Tags::Tag::Reader {
public:
  typedef Tag Reads;

  Reader() = default;
  inline explicit Reader(::capnp::_::StructReader base): _reader(base) {}

  inline ::capnp::MessageSize totalSize() const {
    return _reader.totalSize().asPublic();
  }

#if !CAPNP_LITE
  inline ::kj::StringTree toString() const {
    return ::capnp::_::structString(_reader, *_capnpPrivate::brand());
  }
#endif  // !CAPNP_LITE

  inline  ::int64_t getTime() const;

  inline  ::uint64_t getChannel() const;

private:
  ::capnp::_::StructReader _reader;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::ToDynamic_;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::_::PointerHelpers;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::List;
  friend class ::capnp::MessageBuilder;
  friend class ::capnp::Orphanage;
};

class Tags::Tag::Builder {
public:
  typedef Tag Builds;

  Builder() = delete;  // Deleted to discourage incorrect usage.
                       // You can explicitly initialize to nullptr instead.
  inline Builder(decltype(nullptr)) {}
  inline explicit Builder(::capnp::_::StructBuilder base): _builder(base) {}
  inline operator Reader() const { return Reader(_builder.asReader()); }
  inline Reader asReader() const { return *this; }

  inline ::capnp::MessageSize totalSize() const { return asReader().totalSize(); }
#if !CAPNP_LITE
  inline ::kj::StringTree toString() const { return asReader().toString(); }
#endif  // !CAPNP_LITE

  inline  ::int64_t getTime();
  inline void setTime( ::int64_t value);

  inline  ::uint64_t getChannel();
  inline void setChannel( ::uint64_t value);

private:
  ::capnp::_::StructBuilder _builder;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::ToDynamic_;
  friend class ::capnp::Orphanage;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::_::PointerHelpers;
};

#if !CAPNP_LITE
class Tags::Tag::Pipeline {
public:
  typedef Tag Pipelines;

  inline Pipeline(decltype(nullptr)): _typeless(nullptr) {}
  inline explicit Pipeline(::capnp::AnyPointer::Pipeline&& typeless)
      : _typeless(kj::mv(typeless)) {}

private:
  ::capnp::AnyPointer::Pipeline _typeless;
  friend class ::capnp::PipelineHook;
  template <typename, ::capnp::Kind>
  friend struct ::capnp::ToDynamic_;
};
#endif  // !CAPNP_LITE

// =======================================================================================

inline bool Tags::Reader::hasTags() const {
  return !_reader.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS).isNull();
}
inline bool Tags::Builder::hasTags() {
  return !_builder.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS).isNull();
}
inline  ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>::Reader Tags::Reader::getTags() const {
  return ::capnp::_::PointerHelpers< ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>>::get(_reader.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS));
}
inline  ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>::Builder Tags::Builder::getTags() {
  return ::capnp::_::PointerHelpers< ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>>::get(_builder.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS));
}
inline void Tags::Builder::setTags( ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>::Reader value) {
  ::capnp::_::PointerHelpers< ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>>::set(_builder.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS), value);
}
inline void Tags::Builder::setTags(::kj::ArrayPtr<const  ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>::Reader> value) {
  ::capnp::_::PointerHelpers< ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>>::set(_builder.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS), value);
}
inline  ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>::Builder Tags::Builder::initTags(unsigned int size) {
  return ::capnp::_::PointerHelpers< ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>>::init(_builder.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS), size);
}
inline void Tags::Builder::adoptTags(
    ::capnp::Orphan< ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>>&& value) {
  ::capnp::_::PointerHelpers< ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>>::adopt(_builder.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS), kj::mv(value));
}
inline ::capnp::Orphan< ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>> Tags::Builder::disownTags() {
  return ::capnp::_::PointerHelpers< ::capnp::List< ::capnp::List< ::Tags::Tag,  ::capnp::Kind::STRUCT>,  ::capnp::Kind::LIST>>::disown(_builder.getPointerField(
      ::capnp::bounded<0>() * ::capnp::POINTERS));
}

inline  ::int64_t Tags::Tag::Reader::getTime() const {
  return _reader.getDataField< ::int64_t>(
      ::capnp::bounded<0>() * ::capnp::ELEMENTS);
}

inline  ::int64_t Tags::Tag::Builder::getTime() {
  return _builder.getDataField< ::int64_t>(
      ::capnp::bounded<0>() * ::capnp::ELEMENTS);
}
inline void Tags::Tag::Builder::setTime( ::int64_t value) {
  _builder.setDataField< ::int64_t>(
      ::capnp::bounded<0>() * ::capnp::ELEMENTS, value);
}

inline  ::uint64_t Tags::Tag::Reader::getChannel() const {
  return _reader.getDataField< ::uint64_t>(
      ::capnp::bounded<1>() * ::capnp::ELEMENTS);
}

inline  ::uint64_t Tags::Tag::Builder::getChannel() {
  return _builder.getDataField< ::uint64_t>(
      ::capnp::bounded<1>() * ::capnp::ELEMENTS);
}
inline void Tags::Tag::Builder::setChannel( ::uint64_t value) {
  _builder.setDataField< ::uint64_t>(
      ::capnp::bounded<1>() * ::capnp::ELEMENTS, value);
}


CAPNP_END_HEADER
