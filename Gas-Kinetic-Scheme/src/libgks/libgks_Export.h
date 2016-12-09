
#ifndef libgks_EXPORT_H
#define libgks_EXPORT_H

#ifdef libgks_BUILT_AS_STATIC
#  define libgks_EXPORT
#  define LIBGKS_NO_EXPORT
#else
#  ifndef libgks_EXPORT
#    ifdef libgks_EXPORTS
        /* We are building this library */
#      define libgks_EXPORT __declspec(dllexport)
#    else
        /* We are using this library */
#      define libgks_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef LIBGKS_NO_EXPORT
#    define LIBGKS_NO_EXPORT 
#  endif
#endif

#ifndef LIBGKS_DEPRECATED
#  define LIBGKS_DEPRECATED __declspec(deprecated)
#endif

#ifndef LIBGKS_DEPRECATED_EXPORT
#  define LIBGKS_DEPRECATED_EXPORT libgks_EXPORT LIBGKS_DEPRECATED
#endif

#ifndef LIBGKS_DEPRECATED_NO_EXPORT
#  define LIBGKS_DEPRECATED_NO_EXPORT LIBGKS_NO_EXPORT LIBGKS_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef LIBGKS_NO_DEPRECATED
#    define LIBGKS_NO_DEPRECATED
#  endif
#endif

#endif
