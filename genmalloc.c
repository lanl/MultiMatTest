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

#include <sys/types.h>
#include <sys/queue.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "genmalloc.h"

#ifndef DEBUG
#define DEBUG 0
#endif

SLIST_HEAD(slist_genmalloc_memory_head, genmalloc_memory_entry) genmalloc_memory_head = SLIST_HEAD_INITIALIZER(genmalloc_memory_head);
struct slist_genmalloc_memory_head *genmalloc_memory_headp;
struct genmalloc_memory_entry {
   void *mem_ptr;
   size_t mem_size;
   char *name;
   SLIST_ENTRY(genmalloc_memory_entry) genmalloc_memory_entries;
} *genmalloc_memory_item;

void *genvector_p(const char *name, int inum, size_t elsize, const char *file, const int line)
{
   void *out;
   size_t mem_size;

   mem_size = inum*elsize;
   out      = (void *)malloc((size_t)inum* elsize);
   genmalloc_memory_add(out, name, mem_size);

   return (out);
}

void genvectorfree_p(void *var, const char *file, const int line)
{
   genmalloc_memory_remove(var);
}

void **genmatrix_p(const char *name, int jnum, int inum, size_t elsize, const char *file, const int line)
{
   void **out;
   size_t mem_size;
   const char ptr_name[40];
   sprintf((char *)ptr_name,"%s ptr",name);
  
   mem_size = jnum*sizeof(void *);
   out      = (void **)malloc(mem_size);
   genmalloc_memory_add(out, ptr_name, mem_size);
  
   mem_size = jnum*inum*elsize;
   out[0]    = (void *)malloc((size_t)jnum*(size_t)inum* elsize);
   genmalloc_memory_add(out[0], name, mem_size);
  
   for (int i = 1; i < jnum; i++) {
      out[i] = out[i-1] + inum*elsize;
   }
  
   return (out);
}

void genmatrixfree_p(void **var, const char *file, const int line)
{
   genmalloc_memory_remove(var[0]);
   genmalloc_memory_remove(var);
}

void ***gentrimatrix_p(const char *name, int knum, int jnum, int inum, size_t elsize, const char *file, const int line)
{
   void ***out;
   size_t mem_size;
   const char ptr_name[40];
   sprintf((char *)ptr_name,"%s ptr",name);

   mem_size  = knum*sizeof(void **);
   out       = (void ***)malloc(mem_size);
   genmalloc_memory_add(out, ptr_name, mem_size);

   mem_size  = knum*jnum*sizeof(void *);
   out[0]    = (void **) malloc(mem_size);
   genmalloc_memory_add(out[0], ptr_name, mem_size);

   mem_size  = knum*jnum*inum*elsize;
   out[0][0] = (void *)malloc((size_t)knum*(size_t)jnum*(size_t)inum* elsize);
   genmalloc_memory_add(out[0][0], name, mem_size);

   for (int k = 0; k < knum; k++)
   {
      if (k > 0)
      {
         out[k] = out[k-1] + jnum*sizeof(void *);
         out[k][0] = out[k-1][0] + (jnum*inum);
      }

      for (int j = 1; j < jnum; j++)
      {
         out[k][j] = out[k][j-1] + inum*elsize;
      }
   }
  
   return (out);
}

void gentrimatrixfree_p(void ***var, const char *file, const int line)
{
   genmalloc_memory_remove(var[0][0]);
   genmalloc_memory_remove(var[0]);
   genmalloc_memory_remove(var);
}

void *genmalloc_memory_add_p(void *malloc_mem_ptr, const char *name, size_t size, const char *file, const int line){
   if (SLIST_EMPTY(&genmalloc_memory_head)) SLIST_INIT(&genmalloc_memory_head);

   genmalloc_memory_item = malloc(sizeof(struct genmalloc_memory_entry));
   genmalloc_memory_item->mem_ptr = malloc_mem_ptr;
   genmalloc_memory_item->mem_size = size;
   genmalloc_memory_item->name = strcpy(malloc(strlen(name)+1),name);
   if (DEBUG) printf("GENMALLOC_MEMORY_ADD: DEBUG -- malloc memory pointer is %p\n",malloc_mem_ptr);

   SLIST_INSERT_HEAD(&genmalloc_memory_head, genmalloc_memory_item, genmalloc_memory_entries);

   return(malloc_mem_ptr);
}

void genmalloc_memory_remove_p(void *malloc_mem_ptr, const char *file, const int line){
   SLIST_FOREACH(genmalloc_memory_item, &genmalloc_memory_head, genmalloc_memory_entries){
      if (genmalloc_memory_item->mem_ptr == malloc_mem_ptr) {
         if (DEBUG) printf("GENMALLOC_MEMORY_REMOVE: DEBUG -- freeing malloc memory pointer %p\n",malloc_mem_ptr);
         free(malloc_mem_ptr);
	 free(genmalloc_memory_item->name);
         SLIST_REMOVE(&genmalloc_memory_head, genmalloc_memory_item, genmalloc_memory_entry, genmalloc_memory_entries);
         free(genmalloc_memory_item);
         break;
      }
   }
}

void genmalloc_MB_memory_report_p(const char *file, const int line){
   SLIST_FOREACH(genmalloc_memory_item, &genmalloc_memory_head, genmalloc_memory_entries){
      printf("GENMALLOC REPORT: %20s \t%12lu Bytes, %7lu KB\n",genmalloc_memory_item->name,
             genmalloc_memory_item->mem_size, genmalloc_memory_item->mem_size/1000);
   }
}

void genmalloc_MiB_memory_report_p(const char *file, const int line){
   SLIST_FOREACH(genmalloc_memory_item, &genmalloc_memory_head, genmalloc_memory_entries){
      printf("GENMALLOC REPORT: %20s \t%12lu Bytes, %7lu KiB\n",genmalloc_memory_item->name,
             genmalloc_memory_item->mem_size, genmalloc_memory_item->mem_size/1024);
   }
}

void genmalloc_MiB_memory_total_p(const char *file, const int line){
   size_t Mem_Total = 0;
   SLIST_FOREACH(genmalloc_memory_item, &genmalloc_memory_head, genmalloc_memory_entries){
      Mem_Total += genmalloc_memory_item->mem_size;
   }
   printf("GENMALLOC TOTAL:          Memory Size: \t%12lu Bytes, %7lu KiB, %5lu MiB\n",Mem_Total,Mem_Total/1024,Mem_Total/1024/1024);
}

void genmalloc_MB_memory_total_p(const char *file, const int line){
   size_t Mem_Total = 0;
   SLIST_FOREACH(genmalloc_memory_item, &genmalloc_memory_head, genmalloc_memory_entries){
      Mem_Total += genmalloc_memory_item->mem_size;
   }
   printf("GENMALLOC TOTAL:          Memory Size: \t%12lu Bytes, %7lu KB, %5lu MB\n",Mem_Total,Mem_Total/1000,Mem_Total/1000/1000);
}

void genmem_free_all_p(const char *file, const int line){
   while (!SLIST_EMPTY(&genmalloc_memory_head)) {
      genmalloc_memory_item = SLIST_FIRST(&genmalloc_memory_head);
      if (DEBUG) printf("GENMEM_FREE_ALL: DEBUG -- freeing genmalloc memory %p\n",genmalloc_memory_item->mem_ptr);
      free(genmalloc_memory_item->mem_ptr);
      free(genmalloc_memory_item->name);
      SLIST_REMOVE_HEAD(&genmalloc_memory_head, genmalloc_memory_entries);
      free(genmalloc_memory_item);
   }
}

