/*
 * Copyright (c) 2017, Los Alamos National Security, LLC.
 * All rights Reserved.
 * 
 * This is the code released under LANL Copyright Disclosure C17041/LA-CC-17-041
 * Copyright 2017.  Los Alamos National Security, LLC. This material was produced
 * under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
 * Laboratory (LANL), which is operated by Los Alamos National Security, LLC
 * for the U.S. Department of Energy. See LICENSE file for details.
 *
 * Released under the New BSD License
 *
 * Bob Robey brobey@lanl.gov and Rao Garimella rao@lanl.gov
 */

#ifdef __cplusplus
extern "C"
{
#endif

/* memory routines */
#define genvector(  name, inum, elsize) \
      ( genvector_p(name, inum, elsize, __FILE__, __LINE__) )
#define genvectorfree(  var) \
      ( genvectorfree_p(var, __FILE__, __LINE__) )
#define genmatrix(  name, jnum, inum, elsize) \
      ( genmatrix_p(name, jnum, inum, elsize, __FILE__, __LINE__) )
#define gentrimatrix(  name, knum, jnum, inum, elsize) \
      ( gentrimatrix_p(name, knum, jnum, inum, elsize, __FILE__, __LINE__) )
#define genmatrixfree(  var) \
      ( genmatrixfree_p(var, __FILE__, __LINE__) )
#define gentrimatrixfree(  var) \
      ( gentrimatrixfree_p(var, __FILE__, __LINE__) )

#define genmalloc_memory_add(  malloc_mem_ptr, name, size) \
      ( genmalloc_memory_add_p(malloc_mem_ptr, name, size, __FILE__, __LINE__) )
#define genmalloc_memory_remove(  malloc_mem_ptr) \
      ( genmalloc_memory_remove_p(malloc_mem_ptr, __FILE__, __LINE__) )
#define genmalloc_MiB_memory_report() \
      ( genmalloc_MiB_memory_report_p(__FILE__, __LINE__) )
#define genmalloc_MiB_memory_total() \
      ( genmalloc_MiB_memory_total_p(__FILE__, __LINE__) )
#define genmalloc_MB_memory_report() \
      ( genmalloc_MB_memory_report_p(__FILE__, __LINE__) )
#define genmalloc_MB_memory_total() \
      ( genmalloc_MB_memory_total_p(__FILE__, __LINE__) )
#define genmem_free_all() \
      ( genmem_free_all_p(__FILE__, __LINE__) ) 


void *genvector_p(const char *name, int inum, size_t elsize, const char *file, const int line);
void genvectorfree_p(void *var, const char *file, const int line);
void **genmatrix_p(const char *name, int jnum, int inum, size_t elsize, const char *file, const int line);
void genmatrixfree_p(void **var, const char *file, const int line);
void ***gentrimatrix_p(const char *name, int knum, int jnum, int inum, size_t elsize, const char *file, const int line);
void gentrimatrixfree_p(void ***var, const char *file, const int line);

void *genmalloc_memory_add_p(void *malloc_mem_ptr, const char *name, size_t size, const char *file, const int line);
void genmalloc_memory_remove_p(void *malloc_mem_ptr, const char *file, const int line);
void genmalloc_MiB_memory_report_p(const char *file, const int line);
void genmalloc_MiB_memory_total_p(const char *file, const int line);
void genmalloc_MB_memory_report_p(const char *file, const int line);
void genmalloc_MB_memory_total_p(const char *file, const int line);
void genmem_free_all_p(const char *file, const int line);

#ifdef __cplusplus
}
#endif

