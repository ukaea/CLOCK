/*                                                                                                        
 * 
 *  original LibPng demo code 
 *  "Copyright 2002-2010 Guillaume Cottenceau.                                                              
 *                                                                                                        
 *  This software may be freely redistributed under the terms                                              
 *  of the X11 license."                                                                                    
 *
 *  This version: simple input/output of 16 bit greyscale .png using libpng library released under GPLv3
 *-----------------------------------------------------------------------------------------------------------------------------------
 *    LibPng in the Culham LOop Counting Kit (CLOCK), a library for automated feature detection in irradiated transmission electron micrographs
 *    Copyright (C) 2024  Daniel Mason

 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.

 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.

 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *-----------------------------------------------------------------------------------------------------------------------------------
 *                                                                          
 *                          
 */                                                                                                       
                                                                                                          
#include <stdlib.h>                                                                                       
#include <stdio.h>                                                                                                                                                                        
#include <stdarg.h>                                                                                       
                                                                                                          
#define PNG_DEBUG 3                                                                                       
#include <png.h>                                                                                          
  
float* read_greyscale_png_file(char* file_name,int* width, int* height  )
{
                 
        *width = 0;
        *height = 0;
                                                                                                                                                     
        /* open file and test for it being a png */                                                       
        FILE *fp = fopen(file_name, "rb");                                                                
        if (!fp) {
            printf("Lib_Greyscale::read_greyscale_png_file() error - file io error %s \n",file_name);
            return NULL;                                                                
        }                 
        unsigned char header[8];    // 8 is the maximum size that can be checked                                   
        fread(header, 1, 8, fp);                                                                          
        if (png_sig_cmp(header, 0, 8)) {
            printf("Lib_Greyscale::read_greyscale_png_file() error - file not .png error %s \n",file_name);
            return NULL;                                                                
        }              
        
                                                                                                          
        /* initialize stuff */                                                                            
        png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);                        
                                                                                                          
        if (!png_ptr) {
            printf("Lib_Greyscale::read_greyscale_png_file() error - file not .png error %s \n",file_name);
            return NULL;                                                                
        }                        
                                                                                                          
        png_infop info_ptr = png_create_info_struct(png_ptr);                                                       
        if (!info_ptr) {
            printf("Lib_Greyscale::read_greyscale_png_file() error - read header error %s \n",file_name);
            return NULL;                                                                
        }                                     
                                                                                                          
        if (setjmp(png_jmpbuf(png_ptr))) {
            printf("Lib_Greyscale::read_greyscale_png_file() error - read header error %s \n",file_name);
            return NULL;                                                                
        }                                                  
                                                                                                          
        png_init_io(png_ptr, fp);                                                                         
        png_set_sig_bytes(png_ptr, 8);                                                                    
                                                                                                          
        png_read_info(png_ptr, info_ptr);                                                                 
                                                                                                          
        *width = png_get_image_width(png_ptr, info_ptr);                                                   
        *height = png_get_image_height(png_ptr, info_ptr);                                                 
        png_byte color_type = png_get_color_type(png_ptr, info_ptr);                                               
        png_byte bit_depth = png_get_bit_depth(png_ptr, info_ptr);                                                 
                                                                                                          
                                         
        png_read_update_info(png_ptr, info_ptr);                                                          
                                                                                                                                                        
                                                         
        /* read file */                                                                                   
        if (setjmp(png_jmpbuf(png_ptr))) {
            printf("Lib_Greyscale::read_greyscale_png_file() error - read image error %s \n",file_name);
            return NULL;                                                                
        }                                                             
        png_bytep *  row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * (*height));                  
        
        
        size_t m = sizeof(float) * (*width)*(*height);
        float* img = (float*) malloc(m);
        
        int x,y;
        for (y=0; y<(*height); y++) {                                                                         
             row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr,info_ptr));               
        }                     
                                                                                     
        png_read_image(png_ptr, row_pointers);  
                                                                                                            
        fclose(fp);     
        
        if (color_type == PNG_COLOR_TYPE_RGB) {            
            // printf( "Lib_Greyscale::read_greyscale_png_file info - PNG_COLOR_TYPE_RGB \n");
            for (y=0; y<(*height); y++) {                                                                         
                 png_byte* row = row_pointers[y];                  
                 for (x=0; x<(*width); x++){
                     png_byte* ptr = &(row[x*3]);      
                     png_byte r = ptr[0];
                     png_byte g = ptr[1];
                     png_byte b = ptr[2]; 
                     img[x+(*width)*y] = (r+g+b)/(3*256.0);                
                 }
            }            
        } else if (color_type == PNG_COLOR_TYPE_RGBA ) {          
            // printf( "Lib_Greyscale::read_greyscale_png_file info - PNG_COLOR_TYPE_RGBA \n");   
            for (y=0; y<(*height); y++) {                                                                         
                 png_byte* row = row_pointers[y];                  
                 for (x=0; x<(*width); x++){
                     png_byte* ptr = &(row[x*4]);      
                     png_byte r = ptr[0];
                     png_byte g = ptr[1];
                     png_byte b = ptr[2]; 
                     img[x+(*width)*y] = (r+g+b)/(3*256.0);                
                 }
            }      
        } else if ( (color_type == PNG_COLOR_TYPE_PALETTE)&&(bit_depth == 8) ) {         
            // printf( "Lib_Greyscale::read_greyscale_png_file info - PNG_COLOR_TYPE_PALETTE(8) \n");    
            for (y=0; y<(*height); y++) {                                                                         
                 png_byte* row = row_pointers[y];                  
                 for (x=0; x<(*width); x++){
                     png_byte* ptr = &(row[x]);      
                     png_byte g = ptr[0];                     
                     img[x+(*width)*y] = g/(256.0);             
                 }
            }                        
        } else if (color_type == PNG_COLOR_TYPE_GRAY) {  
            // printf( "Lib_Greyscale::read_greyscale_png_file info - PNG_COLOR_TYPE_GRAY(%d) \n",bit_depth);    
            if (bit_depth<=8){                           
                for (y=0; y<(*height); y++) {                                                                         
                     png_byte* row = row_pointers[y];                  
                     for (x=0; x<(*width); x++){
                         png_byte* ptr = &(row[x]);      
                         png_byte g = ptr[0];                     
                         img[x+(*width)*y] = g/(256.0);                
                     }
                }      
            } else {                
                for (y=0; y<(*height); y++) {                                                                         
                     png_byte* row = row_pointers[y];                  
                     for (x=0; x<(*width); x++){
                         png_byte* ptr = &(row[2*x]);      
                         png_byte g = ptr[0];    
                         png_byte g2 = ptr[1];                     
                         int i = (int) (g*256 + g2);
                         img[x+(*width)*y] = i/65536.0; // g/(256.0) + g2/(16384.0);

                        // if ( (x==0)&&(y==0) ) {
                        //  printf("Lib_Greyscale::read_greyscale_png_file %d,%d,%d,%.8f,%.8f \n",g,g2,i,i/65536.0,img[x+(*width)*y] );
                        // }

                         
                     }
                }   
            }
        } else {
            printf("Lib_Greyscale::read_greyscale_png_file() error -  color_type & bit depth not coded %s \n",file_name);
            return NULL;
        }   
        
        /* cleanup heap allocation */                                                                     
        for (y=0; y<(*height); y++)                                                                          
                free(row_pointers[y]);                                                                    
        free(row_pointers);                                                                               
        
        return img;                                                                                  
}                                                                                                         
                                                 
                                                                                                   
                                                                                                          
void write_greyscale_png_file(char* file_name, int width, int height, float* img )                                                                      
{            
                                                                                                             
        
    
                                                                                                 
        /* create file */                                                                                 
        FILE *fp = fopen(file_name, "wb");                                                                
        if (!fp) {
            printf("Lib_Greyscale::write_greyscale_png_file() error - file io error %s \n",file_name);
            return;
        }   
                                                                                                          
                                                                                                          
        /* initialize stuff */                                                                            
        png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);                       
                                                                                                          
        if (!png_ptr) {
            printf("Lib_Greyscale::write_greyscale_png_file() error - write header error %s \n",file_name);
            return;
        }         
                                                                                                          
        png_infop info_ptr = png_create_info_struct(png_ptr);                                                       
        if (!info_ptr){
            printf("Lib_Greyscale::write_greyscale_png_file() error - write header error %s \n",file_name);
            return;
        }                                        
                                                                                                          
        if (setjmp(png_jmpbuf(png_ptr))) {
            printf("Lib_Greyscale::write_greyscale_png_file() error - write header error %s \n",file_name);
            return;
        }                                            
                                                                                                          
        png_init_io(png_ptr, fp);                                                                         
                                                                                                          
                                                                                                          
        /* write header */                                                                                
        if (setjmp(png_jmpbuf(png_ptr))) {
            printf("Lib_Greyscale::write_greyscale_png_file() error - write header error %s \n",file_name);
            return;
        }                                 
                                                                      
        png_byte bit_depth = (png_byte)16;
        png_byte color_type = PNG_COLOR_TYPE_GRAY;
         
                                            
        png_set_IHDR(png_ptr, info_ptr, width, height,                                                    
                     bit_depth, color_type, PNG_INTERLACE_NONE,                                           
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);                                    
                                                                                                          
        png_write_info(png_ptr, info_ptr);                                                                
                                                                                                          

        /* write bytes */                                                                                 
        if (setjmp(png_jmpbuf(png_ptr))) {
            printf("Lib_Greyscale::write_greyscale_png_file() error - write image error %s \n",file_name);
            return;
        }                                               
        
        
        // convert img data into row_pointers
        png_bytep * row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);             
        int x,y;                      
        for (y=0; y<height; y++) {                                                                         
             row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr,info_ptr));      
        }
        
        
        for (y=0; y<height; y++) {           
             png_byte* row = row_pointers[y];           
             for (x=0; x<width; x++){
                 float f = img[x+width*y];

                 if (f>0.999984741){
                     f = 0.999984741;                   //      1 - 1/65536
                 } else if (f<0.0){
                     f = 0.0;
                 }                                  
                 int i = (int)(f*65536); 
                 png_byte hi = (png_byte)( i>>8 );
                 png_byte lo = (png_byte)( i%256 );
                 
                 png_byte* ptr = &(row[x*2]); 
                 ptr[0] = hi;
                 ptr[1] = lo;


                // if ( (x==0)&&(y==0) ) {
                //     printf("Lib_Greyscale::write_greyscale_png_file %.8f,%d,%d,%d \n",f,i,hi,lo );
                // }


             }
        }
                                                                                                          
        png_write_image(png_ptr, row_pointers);                                                           
                                                                                                          
                                                                                                          
        /* end write */                                                                                   
        if (setjmp(png_jmpbuf(png_ptr))) {
            printf("Lib_Greyscale::write_greyscale_png_file() error - write image error %s \n",file_name);
            return;
        }                                                   
                                                                                                          
        png_write_end(png_ptr, NULL);                                                                     
                                                                                                          
        /* cleanup heap allocation */                                                                     
        for (y=0; y<height; y++)                                                                          
                free(row_pointers[y]);                                                                    
        free(row_pointers);                                                                               
                                                                                                          
        fclose(fp);    
        
        return;                                                                                   
}                                                                                                         
                                            


                
                                                                                                   
                                                                                                          
void write_rgb_png_file(char* file_name, int width, int height, float* img )                                                                      
{            
                                                                                                             
        
    
                                                                                                 
        /* create file */                                                                                 
        FILE *fp = fopen(file_name, "wb");                                                                
        if (!fp) {
            printf("Lib_rgb::write_rgb_png_file() error - file io error %s \n",file_name);
            return;
        }   
                                                                                                          
                                                                                                          
        /* initialize stuff */                                                                             
        png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);                       
                                                                                                          
        if (!png_ptr) {
            printf("Lib_rgb::write_rgb_png_file() error - write header error %s \n",file_name);
            return;
        }         
                                                                                                          
        png_infop info_ptr = png_create_info_struct(png_ptr);                                                       
        if (!info_ptr){
            printf("Lib_rgb::write_rgb_png_file() error - write header error %s \n",file_name);
            return;
        }                                        
                                                                                                          
        if (setjmp(png_jmpbuf(png_ptr))) {
            printf("Lib_rgb::write_rgb_png_file() error - write header error %s \n",file_name);
            return;
        }                                            
                                                                                                          
        png_init_io(png_ptr, fp);                                                                         
                                                                                                          
                                                                                                          
        /* write header */                                                                                
        if (setjmp(png_jmpbuf(png_ptr))) {
            printf("Lib_rgb::write_rgb_png_file() error - write header error %s \n",file_name);
            return;
        }                                 
                                                                      
        png_byte bit_depth = (png_byte)8;
        png_byte color_type = PNG_COLOR_TYPE_RGB;
         
                                            
        png_set_IHDR(png_ptr, info_ptr, width, height,                                                    
                     bit_depth, color_type, PNG_INTERLACE_NONE,                                           
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);                                     
                                                                                                          
        png_write_info(png_ptr, info_ptr);                                                                
                                                                                                          
  
        /* write bytes */                                                                                 
        if (setjmp(png_jmpbuf(png_ptr))) {
            printf("Lib_rgb::write_rgb_png_file() error - write image error %s \n",file_name);
            return;
        }                                               
        
        
        // convert img data into row_pointers   
        png_bytep * row_pointers = (png_bytep*) malloc( sizeof(png_bytep) * height);             
        int x,y;                      
        for (y=0; y<height; y++) {                                                                         
             row_pointers[y] = (png_byte*) malloc( png_get_rowbytes(png_ptr,info_ptr) );      
        }
        
        
        for (y=0; y<height; y++) {           
             png_byte* row = row_pointers[y];           
             for (x=0; x<width; x++){
                 
                 float r = img[ 3*x+0 + 3*width*y ];
                 float g = img[ 3*x+1 + 3*width*y ];
                 float b = img[ 3*x+2 + 3*width*y ];

                 r = (r>0.99609375)?0.99609375:(r<0?0:r);           //  255/256
                 g = (g>0.99609375)?0.99609375:(g<0?0:g);
                 b = (b>0.99609375)?0.99609375:(b<0?0:b);
                 
                                                        
                 
                 png_byte* ptr = &(row[x*3]); 
                 ptr[0] = (png_byte)( (int)(r*256) );
                 ptr[1] = (png_byte)( (int)(g*256) );
                 ptr[2] = (png_byte)( (int)(b*256) );
                  
             }
        }
                                                                                                          
        png_write_image(png_ptr, row_pointers);                                                           
                                                                                                          
                                                                                                          
        /* end write */                                                                                   
        if (setjmp(png_jmpbuf(png_ptr))) {
            printf("Lib_rgb::write_rgb_png_file() error - write image error %s \n",file_name);
            return;
        }                                                   
                                                                                                          
        png_write_end(png_ptr, NULL);                                                                     
                                                                                                          
        /* cleanup heap allocation */                                                                     
        for (y=0; y<height; y++)                                                                          
                free(row_pointers[y]);                                                                    
        free(row_pointers);                                                                               
                                                                                                          
        fclose(fp);    
        
        return;                                                                                   
}                                                                                                          

void destroy_png_storage(float *img)
{
   free(img);
}