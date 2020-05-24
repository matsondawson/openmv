/*
 * This file is part of the OpenMV project.
 *
 * By Matt Dawson
 * matsondawson@gmail.com
 * 
 * References:
 * 1/ https://docs.microsoft.com/en-us/typography/opentype/spec/kern
 *
 * This work is licensed under the MIT license, see the file LICENSE for details.
 *
 */
#include "imlib.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "py_assert.h"

#define ttf_error(msg) PY_ASSERT_TRUE_MSG(true, msg);
#define ttf_assert(cond,msg) PY_ASSERT_TRUE_MSG(cond,msg)

#define ttf_fixed int32_t
#define ttf_fword int16_t
#define ttf_ufword uint16_t

const int ON_CURVE = 1, X_IS_BYTE = 2, Y_IS_BYTE = 4, REPEAT = 8, X_DELTA = 16, Y_DELTA = 32;

typedef struct {
    // A tag to indicate the OFA scaler to be used to rasterize this font; see the note on the scaler type below for more information.
    uint32_t	scaler_type; // Not used
    // number of tables
    uint16_t	numTables;
    // (maximum power of 2 <= numTables)*16
    uint16_t	searchRange; // Not used
    // log2(maximum power of 2 <= numTables)
    uint16_t	entrySelector; // Not used
    // numTables*16-searchRange
    uint16_t	rangeShift; // Not used
} ttf_offset_subtable_t;

typedef struct {
    // 4-byte identifier
    uint32_t tag;
    // checksum for this table
    uint32_t	checkSum;
    // offset from beginning of sfnt
    uint32_t	offset;
    // length of this table in byte (actual length not padded length)
    uint32_t	length;
} ttf_table_directory;

typedef struct {
    // 0x00010000 if (version 1.0)
    ttf_fixed	version;
    // set by font manufacturer
    ttf_fixed	fontRevision;
    // To compute: set it to 0, calculate the checksum for the 'head' table and put it in the table directory, sum the entire font as uint32, then store B1B0AFBA - sum. The checksum for the 'head' table will not be wrong. That is OK.
    uint32_t	checkSumAdjustment;
    //	set to 0x5F0F3CF5
    uint32_t	magicNumber;
    uint16_t	flags;
    /*
    bit 0 - y value of 0 specifies baseline
    bit 1 - x position of left most black bit is LSB
    bit 2 - scaled point size and actual point size will differ (i.e. 24 point glyph differs from 12 point glyph scaled by factor of 2)
    bit 3 - use integer scaling instead of fractional
    bit 4 - (used by the Microsoft implementation of the TrueType scaler)
    bit 5 - This bit should be set in fonts that are intended to e laid out vertically, and in which the glyphs have been drawn such that an x-coordinate of 0 corresponds to the desired vertical baseline.
    bit 6 - This bit must be set to zero.
    bit 7 - This bit should be set if the font requires layout for correct linguistic rendering (e.g. Arabic fonts).
    bit 8 - This bit should be set for a GX font which has one or more metamorphosis effects designated as happening by default.
    bit 9 - This bit should be set if the font contains any strong right-to-left glyphs.
    bit 10 - This bit should be set if the font contains Indic-style rearrangement effects.
    bits 11-12 - Defined by Adobe.*/
    // range from 64 to 16384
    uint16_t	unitsPerEm;
    uint32_t created_low, created_high;
    uint32_t modified_low, modified_high;
    //ttf_long_date_time created; // created	international date
    //ttf_long_date_time modified; // modified	international date
    ttf_fword	xMin; // for all glyph bounding boxes
    ttf_fword	yMin; // for all glyph bounding boxes
    ttf_fword	xMax; // for all glyph bounding boxes
    ttf_fword	yMax; // for all glyph bounding boxes
    uint16_t	macStyle;
    /*
    bit 0 bold
    bit 1 italic
    bit 2 underline
    bit 3 outline
    bit 4 shadow
    bit 5 condensed (narrow)
    bit 6 extended
    */
    // smallest readable size in pixels
    uint16_t	lowestRecPPEM;
    int16_t	fontDirectionHint;
    /*
    0 Mixed directional glyphs
    1 Only strongly left to right glyphs
    2 Like 1 but also contains neutrals
    -1 Only strongly right to left glyphs
    -2 Like -1 but also contains neutrals*/
    int16_t	indexToLocFormat;//	0 for short offsets, 1 for long
    int16_t	glyphDataFormat;//	0 for current format
} ttf_table_head;

typedef struct {
    // If the number of contours is positive or zero, it is a single glyph;
    // If the number of contours is -1, the glyph is compound
    int16_t numberOfContours;
    // Minimum x for coordinate data
    ttf_fword	xMin;
    // Minimum y for coordinate data	
    ttf_fword	yMin;
    // Maximum x for coordinate data
    ttf_fword	xMax;
    // Maximum y for coordinate data
    ttf_fword	yMax;
} ttf_table_glyf;

typedef struct {
    // Set to 0
    uint16_t format;
    // Length in bytes of the subtable (set to 262 for format 0)
    uint16_t length;
    // Language code for this encoding subtable, or zero if language-independent
    uint16_t language;
    // An array that maps character codes to glyph index values
    uint8_t glyphIndexArray[256];
} table_cmap_format_0;

typedef struct {
    uint16_t	format; //	Format number is set to 4	 
    uint16_t	length; //	Length of subtable in bytes	 
    uint16_t	language; //	Language code for this encoding subtable, or zero if language-independent	 
    uint16_t	segCountX2; //	2 * segCount	 
    uint16_t	searchRange; //	2 * (2**FLOOR(log2(segCount)))	 
    uint16_t	entrySelector; //	log2(searchRange/2)
    uint16_t	rangeShift; //	(2 * segCount) - searchRange	 
    uint16_t	endCode0; // * segCount; //	Ending character code for each segment, last = 0xFFFF.
    /*uint16_t	reservedPad; //	This value should be zero	
    uint16_t	startCode[segCount]	; //Starting character code for each segment	
    uint16_t	idDelta[segCount]; //	Delta for all character codes in segment	 
    uint16_t	idRangeOffset[segCount]; //	Offset in bytes to glyph indexArray, or 0	 
    uint16_t	glyphIndexArray[variable]; //	Glyph index array*/
} table_cmap_format_4;

typedef struct {
    // Platform identifier
    uint16_t platformID;
    // Platform-specific encoding identifier
    uint16_t platformSpecificID;
    // Offset of the mapping table
    uint32_t offset;
} ttf_table_cmap_subtable;

typedef struct {
    uint16_t version;
    uint16_t numberSubtables;
    ttf_table_cmap_subtable subtable[];
} ttf_table_cmap;

typedef struct {
    // 0x00010000 (1.0)
    ttf_fixed	version;
    // Distance from baseline of highest ascender
    ttf_fword	ascent;
    // Distance from baseline of lowest descender
    ttf_fword	descent;
    // typographic line gap
    ttf_fword	lineGap;
    // must be consistent with horizontal metrics
    ttf_ufword advanceWidthMax;
    // must be consistent with horizontal metrics
    ttf_fword	minLeftSideBearing;
    // must be consistent with horizontal metrics
    ttf_fword	minRightSideBearing;
    // max(lsb + (xMax-xMin))
    ttf_fword	xMaxExtent;
    // used to calculate the slope of the caret (rise/run) set to 1 for vertical caret
    int16_t	caretSlopeRise;
    // 0 for vertical
    int16_t	caretSlopeRun;
    // set value to 0 for non-slanted fonts
    ttf_fword	caretOffset;
    int16_t	reserved[4];
    // 0 for current format
    int16_t	metricDataFormat;
    // number of advance widths in metrics table
    uint16_t numOfLongHorMetrics;
} ttf_table_hhea;

typedef struct {
    uint16_t advanceWidth;
    int16_t leftSideBearing;
} longHorMetric;

typedef struct {
    const ttf_table_head* head;
    const void* loca;
    const ttf_table_glyf *glyf;
    const ttf_table_cmap *cmap;
    const ttf_table_hhea *hhea;
    const longHorMetric *longHorMetric;

    union {
        const void *data;
        const ttf_offset_subtable_t* offset_subtable;
    };
} ttf_font_instance;


typedef struct ttf_point {
    float x, y;
    uint8_t flags;
} ttf_point;

typedef struct {
    uint16_t numPoints;
    uint16_t *endPtsOfContours;//	Array of last points of each contour; n is the number of contours; array entries are point indices
    uint16_t instructionLength;//	Total number of bytes needed for instructions
    uint8_t	*instructions; // Array of instructions for this glyph
    ttf_point	*coordinates;
} ttf_glyph_data;


typedef struct {
    float x;
    int8_t d;
} ttf_intersect_list_element;

typedef struct {
    uint32_t index;
    uint32_t count;
    ttf_intersect_list_element list[32];
} ttf_intersect_list;


/************************************************* truetypefont object functions *************************************************/

const void* ttf_find_table(ttf_font_instance* ttf, char *tag) {
    ttf_table_directory* table_directory = (ttf_table_directory*)(ttf->data + sizeof(ttf_offset_subtable_t));
    for(int i = 0; i < __REV16(ttf->offset_subtable->numTables); i++) {
        if (table_directory->tag == *(uint32_t*)tag) {
            // printf("table %s offset %lu\n", tag, __REV(table_directory->offset));
            return ttf->data + __REV(table_directory->offset);
        }
        table_directory++;
    }
    return NULL;
}

void ttf_init(ttf_font_instance* ttf, const uint8_t *data) {
    ttf->data = data;
    ttf->head = ttf_find_table(ttf, "head");
    ttf_assert(ttf->head != NULL, "No head table");
    ttf->loca = ttf_find_table(ttf, "loca");
    ttf_assert(ttf->loca != NULL, "No loca table");
    ttf->glyf = ttf_find_table(ttf, "glyf");
    ttf_assert(ttf->glyf != NULL, "No glyf table");
    ttf->cmap = ttf_find_table(ttf, "cmap");
    ttf_assert(ttf->cmap != NULL, "No cmap table");
    ttf->hhea = ttf_find_table(ttf, "hhea");
    ttf_assert(ttf->hhea != NULL, "No hhea table");
    ttf->longHorMetric = ttf_find_table(ttf, "hmtx");
    ttf_assert(ttf->longHorMetric != NULL, "No hmtx table");
}

longHorMetric ttf_get_horizontal_metrics(ttf_font_instance* ttf, uint32_t glyphIndex) {
    // If glyph exists in table
    uint16_t numOfLongHorMetrics = __REV16(ttf->hhea->numOfLongHorMetrics);
    
    return ttf->longHorMetric[(glyphIndex < numOfLongHorMetrics) ? glyphIndex : (numOfLongHorMetrics - 1)];
}

/**
 * Find glyph in the glyph table using hte loca table.
 *
 * @param ttf Truetype font instance.
 * @param glyph_index Index of glyph.
 * @return Pointer to the glyph in the 'glyf' table if found else NULL
 */
const ttf_table_glyf *ttf_get_glyph_ptr(ttf_font_instance* ttf, uint32_t glyph_index) {
    uint32_t offset, next;

    if (__REV16(ttf->head->indexToLocFormat) == 1) {
        uint32_t *ptr = ((uint32_t*)(ttf->loca)) + glyph_index;
        offset = __REV(*(ptr++));
        next = __REV(*ptr);
    } else {
        uint16_t *ptr = (uint16_t*)(ttf->loca) + glyph_index;
        offset = __REV16(*(ptr++)) << 1;
        next = __REV16(*ptr) << 1;
    }

    // indicates glyph has no outline( eg space)
    if (offset == next) return NULL;

    return (void*)ttf->glyf + offset;
}

/**
 * Read an interpret the glyph from the 'glyf' table.
 *
 * @param ttf Truetype font instance.
 * @param glyph_data Where to put the interpreted data.
 * @param glyph_index Index of glyph.
 */
bool ttf_read_glyph(ttf_font_instance* ttf, ttf_glyph_data *glyph_data, uint32_t glyph_index) {
    const ttf_table_glyf *glyph = ttf_get_glyph_ptr(ttf, glyph_index);

    // Glyph is not in the font table.
    if (!glyph) return false;

    // Glyph is probably space
    if (!glyph->numberOfContours) return false;

    glyph_data->endPtsOfContours = (void*)glyph + sizeof(ttf_table_glyf);
    
    uint16_t numberOfContours = __REV16(glyph->numberOfContours);
    
    // Total number of points is the max of endPtsOfContours+1
    int numPoints = -1;
    for(int i = 0; i < numberOfContours; i++) {
        uint16_t v = __REV16(glyph_data->endPtsOfContours[i]);
        if (v > numPoints) numPoints = v;
    }
    numPoints++;

    // Allocate space for interpreted coordinates, this is freed in the function calling this
    glyph_data->coordinates = fb_alloc(numPoints * sizeof(ttf_point), FB_ALLOC_NO_HINT);
    glyph_data->numPoints = numPoints;
    glyph_data->instructionLength = __REV16( *(glyph_data->endPtsOfContours + numberOfContours) );
    glyph_data->instructions = (uint8_t*)(glyph_data->endPtsOfContours + numberOfContours + 1);

    uint8_t *flags = (glyph_data->instructions + glyph_data->instructionLength);

    if (numberOfContours == -1) {
        ttf_error("Compound glyphs not supported");
    } else {
        // Copy flags to coordinates
        int i=0;
        for (int point = 0; point < numPoints;) {
            uint8_t flag = flags[i++];
            glyph_data->coordinates[point].x = 0;
            glyph_data->coordinates[point].y = 0;
            glyph_data->coordinates[point++].flags  = flag;

            if (flag & REPEAT) {
                int repeatCount = flags[i++];
                while (repeatCount--) {
                    glyph_data->coordinates[point].x = 0;
                    glyph_data->coordinates[point].y = 0;
                    glyph_data->coordinates[point++].flags  = flag;
                }
            }
        }
        
        uint8_t *encoded_coordinate_data = ((void*)flags) + i;

        // Decode x coordinates
        {
            int16_t value = 0;

            for (int i = 0; i < numPoints; i++) {
                uint8_t flag = glyph_data->coordinates[i].flags;
                if (flag & X_IS_BYTE) {
                    uint8_t d = *(encoded_coordinate_data++);
                    if (flag & X_DELTA) value += d; else  value -= d;
                } else if (~flag & X_DELTA) {
                    value += __REV16(*(int16_t*)encoded_coordinate_data);
                    encoded_coordinate_data += 2;
                }

                glyph_data->coordinates[i].x = value;
            }
        }

        // Decode y coordinates
        {
            int16_t value = 0;

            for (int i = 0; i < numPoints; i++) {
                uint8_t flag = glyph_data->coordinates[i].flags;
                if (flag & Y_IS_BYTE) {
                    uint8_t d = *(encoded_coordinate_data++);
                    if (flag & Y_DELTA) value += d; else  value -= d;
                } else if (~flag & Y_DELTA) {
                    value += __REV16(*(int16_t*)encoded_coordinate_data);
                    encoded_coordinate_data += 2;
                }

                glyph_data->coordinates[i].y = value;
            }
        }
    }

    return true;
}

/**
 * Convert character code to glyph index.
 *
 * @param ttf Truetype font instance.
 * @param char_code Character code.
 * @return Index of glyph.
 */
int ttf_table_cmap_char_to_index(ttf_font_instance* ttf, uint16_t char_code) {
    for (int i = 0; i < __REV16(ttf->cmap->numberSubtables); i++) {
        // platforms are: 
        // 0 - Unicode -- use specific id 6 for full coverage. 0/4 common.
        // 1 - MAcintosh (Discouraged)
        // 2 - reserved
        // 3 - Microsoft
        const ttf_table_cmap_subtable *subtable = &(ttf->cmap->subtable[i]);
        uint16_t platformID = __REV16(subtable->platformID);
        uint16_t platformSpecificID = __REV16(subtable->platformSpecificID);
        uint32_t offset = __REV(subtable->offset);

        if (platformID == 3 && (platformSpecificID <= 1)) {
            void *cmap_table_ptr = (void*)(ttf->cmap) + offset;
            uint16_t format = __REV16(*(uint32_t*)cmap_table_ptr);

            if (format == 0) {
                table_cmap_format_0 *cmap0 = cmap_table_ptr;
                if (char_code >= 0 && char_code < 256) {
                    return cmap0->glyphIndexArray[char_code];
                }
            } else if (format == 4) {
                table_cmap_format_4 *cmap4 = cmap_table_ptr;

                uint16_t segCount = __REV16(cmap4->segCountX2) >> 1;
                uint16_t *endCode_ptr= &(cmap4->endCode0);
                uint16_t *startCode_ptr= endCode_ptr + segCount + 1;
                uint16_t *idDelta_ptr= startCode_ptr + segCount;
                uint16_t *idRangeOffset_ptr= idDelta_ptr + segCount;
                // uint16_t *glyphIndexArray_ptr = idRangeOffset_ptr + segCount;
                
                for(int segment_index=0; segment_index < segCount; segment_index++) {
                    if (__REV16(startCode_ptr[segment_index]) <= char_code && __REV16(endCode_ptr[segment_index]) >= char_code) {
                        uint16_t index;
                        if (idRangeOffset_ptr[segment_index]) {
                            uint16_t *glyphIndexAddress = &(idRangeOffset_ptr[segment_index]) + (__REV16(idRangeOffset_ptr[segment_index]) / 2) + (char_code - __REV16(startCode_ptr[segment_index]));
                            index = __REV16(*glyphIndexAddress);
                        } else {
                            index = (__REV16(idDelta_ptr[segment_index]) + char_code) & 0xffff;
                        }

                        return index;
                    }
                }
            }
            else {
                ttf_error("Unsupported cmap format");
            }
        }
        else {
            ttf_error("Unsupported cmap platformId");
        }
    }
    
    return 0;
}

void ttf_intersect_list_init(ttf_intersect_list *intersect_list) {
    intersect_list->index = intersect_list->count = 0;
}

void ttf_intersect_list_add(ttf_intersect_list *intersect_list, float x, bool d) {
    if (intersect_list->count < sizeof(intersect_list->list) / sizeof(ttf_intersect_list_element)) {
        ttf_intersect_list_element *element = &intersect_list->list[intersect_list->count++];
        element->x = x;
        element->d = d ? 1 : -1;
    }
}

void ttf_test_and_add_intersection(ttf_intersect_list *intersect_list, int32_t y, ttf_point *point1, ttf_point *point2) {
    float y1 = point1->y, y2 = point2->y;
    bool possible = (y1 != y2) && ((y >= y1 && y < y2) || (y < y1 && y >= y2));

    if (possible) {
        float x1 = point1->x, x2 = point2->x;
        float g = (x1 - x2) / (y1 - y2);
        float x = g * (y - y1) + x1;
        ttf_intersect_list_add(intersect_list, x, y1 > y2);
    }
    else {
        // This adds horizontal lines that are missed by scan line.
        if (y1 == y2 && floorf(y1) == floorf(y2) && floorf(y1) == y) {
            float x1 = point1->x, x2 = point2->x;
            if (x1 < x2) {
                ttf_intersect_list_add(intersect_list, x1, false);
                ttf_intersect_list_add(intersect_list, x2, true);
            }
        }
    }
}

void ttf_intersect_list_sort(ttf_intersect_list *intersect_list) {
    // Don't look at this bubble sort :)
    // There's usually no more than 4 values so it doesn't really matter.
    ttf_intersect_list_element *list = intersect_list->list;
    int32_t count = intersect_list->count;

    for(int i=0; i < count - 1; i++) {
        for(int j=0; j < count - i - 1; j++) {
            if (list[j].x > list[j + 1].x || (list[j].x == list[j + 1].x && list[j].d > list[j + 1].d)) {
                float x = list[j].x;
                list[j].x = list[j + 1].x;
                list[j + 1].x = x;

                int8_t d = list[j].d;
                list[j].d = list[j + 1].d;
                list[j + 1].d = d;
            }
        }
    }
}

ttf_intersect_list_element *ttf_intersect_list_get(ttf_intersect_list *intersect_list) {
    if (intersect_list->index < intersect_list->count) {
        return &intersect_list->list[intersect_list->index];
    }

    return NULL;
}

bool ttf_intersect_list_next(ttf_intersect_list *intersect_list) {
    if (intersect_list->index < intersect_list->count) {
        intersect_list->index++;
        if (intersect_list->index < intersect_list->count) return true;
    }
    
    return false;
}

void ttf_grid_fit(ttf_font_instance* ttf, ttf_glyph_data *glyph_data) {
    float longest_x = 10000, max_len = 0;
    float min_x = 10000, max_x = -10000;

    for (int i = 0; i < glyph_data->numPoints; i++) {
        float x2 = glyph_data->coordinates[i].x;

        min_x = IM_MIN(x2, min_x);
        max_x = IM_MAX(x2, max_x);

        if (i > 0) {
            float x1 = glyph_data->coordinates[i-1].x;

            if (x1 == x2) {
                float y1 = glyph_data->coordinates[i-1].y;
                float y2 = glyph_data->coordinates[i].y;

                float len = fabs(y1-y2);
                if (len > max_len) {
                    longest_x = x1;
                    max_len = len;
                }
            }
        }
    }

    float shift_x = longest_x - floorf(longest_x);

    for (int i = 0; i < glyph_data->numPoints; i++) {
        glyph_data->coordinates[i].x -= shift_x;
    }
    longest_x = floorf(longest_x);

    float width = max_x - min_x;
    float x_mid = min_x + (width / 2);
    if (longest_x > x_mid) {
        float lhs_width = longest_x - min_x;
        float lhs_width_ceil = longest_x - ceilf(min_x);

        // if longest_x was on rhs, then shrink the lhs most point to pixel
        float scale = lhs_width_ceil / lhs_width;
        
        for (int i = 0; i < glyph_data->numPoints; i++) {
            glyph_data->coordinates[i].x = longest_x - (longest_x - glyph_data->coordinates[i].x) * scale;
            // lx = 10, x= 5
            // 10 - (10 - 5) = 5
        }
    } else {
        // if min_x was on rhs, then expand lhs most point to pixel
        float rhs_width = max_x - longest_x;
        float rhs_width_floor = floorf(max_x) - longest_x;

        // if longest_x was on rhs, then shrink the lhs most point to pixel
        float scale = rhs_width_floor / rhs_width;
        
        for (int i = 0; i < glyph_data->numPoints; i++) {
            glyph_data->coordinates[i].x = longest_x + (glyph_data->coordinates[i].x - longest_x) * scale;
        }
    }
}

ttf_point last_draw_point;

void ttf_move_to(ttf_point *point) {
    last_draw_point = *point;
}

void ttf_line_to(ttf_intersect_list *intersect_list, float y, ttf_point *point) {
    ttf_test_and_add_intersection(intersect_list, y, &last_draw_point, point);
    last_draw_point = *point;
}

void ttf_curve_to(ttf_intersect_list *intersect_list, float y, ttf_point *control_point, ttf_point *point2) {
    float y1 = last_draw_point.y, y2 = control_point->y;
    bool possible = (floorf(y1) == floorf(y2) && floorf(y) == floorf(y2)) || ((y >= y1 && y < y2) || (y < y1 && y >= y2));
    if (!possible) {
        y1 = control_point->y, y2 = point2->y;
        possible = (floorf(y1) == floorf(y2) && floorf(y) == floorf(y2)) || ((y >= y1 && y < y2) || (y < y1 && y >= y2));
    }
    if (!possible) {
        last_draw_point = *point2;
        return;
    }

    float max_distance = IM_MAX(fabs(last_draw_point.x - point2->x), fabs(last_draw_point.y - point2->y));
    // At least one curve point every 4 pixels
    int steps = ceilf(max_distance / 2);
    steps = (steps + 1) & ~1;
    
    if (steps <= 2) {
        ttf_test_and_add_intersection(intersect_list, y, &last_draw_point, control_point);
        ttf_test_and_add_intersection(intersect_list, y, control_point, point2);
        last_draw_point = *point2;
    } else {
        ttf_point point0 = last_draw_point;
        ttf_point temp;
        for(int i = 0; i <= steps; i++) {
            float t = (float)i / steps;
            float t2 = t * t;
            float t1m2 = (1 - t);
            t1m2 *= t1m2;

            temp.x = control_point->x + t1m2 * (point0.x - control_point->x) + t2 * (point2->x - control_point->x); 
            temp.y = control_point->y + t1m2 * (point0.y - control_point->y) + t2 * (point2->y - control_point->y);

            ttf_test_and_add_intersection(intersect_list, y, &last_draw_point, &temp);

            last_draw_point = temp;
        }
    }
}

float ttf_draw_glyph(ttf_font_instance* ttf, image_t *img, uint16_t glyph_index, uint32_t color, int32_t location_x, int32_t location_y, float size_pixels, int32_t offset_x_units, int32_t offset_y_units) {
    const ttf_table_glyf *glyf = ttf_get_glyph_ptr(ttf, glyph_index);
    
    // glyph not found, probably space character
    if (!glyf) return 0;

    // Clean and convert color to pixel format
    uint32_t fr = 0, fg = 0, fb = 0;

    switch(img->bpp) {
        case IMAGE_BPP_BINARY:
        case IMAGE_BPP_GRAYSCALE:
            color &= 255;
            break;

        case IMAGE_BPP_RGB565:
            color &= 65535;
            fr = COLOR_RGB565_TO_R5(color);
            fg = COLOR_RGB565_TO_G6(color);
            fb = COLOR_RGB565_TO_B5(color);
            break;
    }

    uint16_t units_per_em = __REV16(ttf->head->unitsPerEm);

    // Y-coordinate is inverted so we draw from the bottom up, hence we move offset down one line
    offset_y_units += units_per_em;

    float pixels_per_unit = size_pixels / units_per_em;

    int32_t offset_x_pixels = floorf(location_x + offset_x_units * pixels_per_unit);
    int32_t offset_y_pixels = floorf(location_y + offset_y_units * pixels_per_unit);
    
    ttf_glyph_data glyph_data;
    if (!ttf_read_glyph(ttf, &glyph_data, glyph_index)) {
        return 0;
    }
    
    // Round?
    int32_t xMin = floorf(((int16_t)__REV16(glyf->xMin)) * pixels_per_unit);
    int32_t xMax = ceilf(((int16_t)__REV16(glyf->xMax)) * pixels_per_unit);
    int32_t yMin = floorf(((int16_t)__REV16(glyf->yMin)) * pixels_per_unit);
    int32_t yMax = ceilf(((int16_t)__REV16(glyf->yMax)) * pixels_per_unit);
    
    if (yMax > offset_y_pixels) {
        yMax = offset_y_pixels;
    }

    if (offset_y_pixels - yMin >= img->h) {
        yMin = offset_y_pixels - img->h + 1;
    }

    if ((offset_x_pixels + xMin) < 0) {
        xMin = -offset_x_pixels;
    }

    if (offset_x_pixels + xMax >= img->w) {
        xMax = (img->w - offset_x_pixels- 1);
    }

    if (xMin > xMax || yMin > yMax) return 0;
    
    for (int i = 0; i < glyph_data.numPoints; i++) {
        glyph_data.coordinates[i].x *= pixels_per_unit;
        glyph_data.coordinates[i].y *= pixels_per_unit;
    }

    ttf_grid_fit(ttf, &glyph_data);

    ttf_intersect_list intersect_list;
    for (int y = yMin; y <= yMax ; y++) {
        ttf_intersect_list_init(&intersect_list);
        
        for (int p = 0, s = 0, c = 0, contour_start = 0; p < glyph_data.numPoints; p++) {
            ttf_point *point = &(glyph_data.coordinates[p]);
            
            if (s == 0) {
                s = 1;
                ttf_move_to(point);
            } else if (s == 1) {
                if (point->flags & ON_CURVE) {
                    ttf_line_to(&intersect_list, y, point);
                }
                else {
                    s = 2;
                }
            } else {
                ttf_point *prev = &(glyph_data.coordinates[p - 1]);
                if (point->flags & ON_CURVE) {
                    ttf_curve_to(&intersect_list, y, prev, point);
                    s = 1;
                }
                else {
                    ttf_point end_point;
                    end_point.x = (prev->x + point->x) / 2;
                    end_point.y = (prev->y + point->y) / 2;
                    ttf_curve_to(&intersect_list, y, prev, &end_point);
                }
            }

            if (p == __REV16(glyph_data.endPtsOfContours[c])) {
                ttf_point *prev = point;

                point = &(glyph_data.coordinates[contour_start]);

                if (point->flags & ON_CURVE) {
                    ttf_curve_to(&intersect_list, y, prev, point);
                } else {
                    ttf_point end_point;
                    end_point.x = (prev->x + point->x) / 2;
                    end_point.y = (prev->y + point->y) / 2;
                    ttf_curve_to(&intersect_list, y, prev, &end_point);
                }

                contour_start = p + 1;
                s = 0;
                c++;
            }
        }
        
        ttf_intersect_list_sort(&intersect_list);
        
        int inside = 0;

        // Remove leading pixels outside of range, but still calculate wheter line is inside the curve.
        {
            ttf_intersect_list_element *intersect = NULL;
            
            while ((intersect = ttf_intersect_list_get(&intersect_list)) && intersect->x < xMin) {
                inside += intersect->d;
                ttf_intersect_list_next(&intersect_list);
            }
        }

        int32_t y_pixel = offset_y_pixels - y;
        const void *row_ptr = imlib_compute_row_ptr(img, y_pixel);
        
        for (int x = xMin; x <= xMax; x++) {
            float last_x = -100000;
            float delta_sum = 0;
            ttf_intersect_list_element *intersect = NULL;
            float this_x;

            while((intersect = ttf_intersect_list_get(&intersect_list)) && x == floorf(this_x = intersect->x)) {
                if (inside < 0) {
                    delta_sum += (this_x - ((last_x == -100000) ? x : last_x));
                }
                inside += intersect->d;
                last_x = this_x;

                if (!ttf_intersect_list_next(&intersect_list)) break;
            }

            if (last_x == -100000) {
                // If there were no x intercepts in this pixel then fill the pixel.
                delta_sum = (inside < 0) ? 1 : 0;
            } else if (inside < 0) {
                // Otherwise fill the remainder of the pixel
                delta_sum += (x + 1.0f - last_x);
            }

            if (delta_sum != 0) {
                int32_t x_pixel = offset_x_pixels + x;
                uint32_t alpha = delta_sum * 255;

                switch(img->bpp) {
                    case IMAGE_BPP_BINARY: {
                        if (alpha >> 7) {
                            IMAGE_PUT_BINARY_PIXEL_FAST((uint32_t*)row_ptr, x_pixel, 1);
                        }
                        break;
                    }
                    case IMAGE_BPP_GRAYSCALE: {
                        if (alpha == 255) {
                            IMAGE_PUT_GRAYSCALE_PIXEL_FAST((uint8_t*)row_ptr, x_pixel, 255 );
                        } else {
                            uint32_t img_pixel = IMAGE_GET_GRAYSCALE_PIXEL_FAST((uint8_t*)row_ptr, x_pixel);

                            IMAGE_PUT_GRAYSCALE_PIXEL_FAST((uint8_t*)row_ptr, x_pixel, ((img_pixel * (256 - alpha))>>8) + alpha);
                        }
                        break;
                    }
                    case IMAGE_BPP_RGB565: {
                        if (alpha==255) {
                            IMAGE_PUT_RGB565_PIXEL_FAST((uint16_t*)row_ptr, x_pixel, color);
                        } else {
                            uint32_t img_pixel = IMAGE_GET_RGB565_PIXEL_FAST((uint16_t*)row_ptr, x_pixel);
                            uint32_t vr = COLOR_RGB565_TO_R5(img_pixel);
                            uint32_t vg = COLOR_RGB565_TO_G6(img_pixel);
                            uint32_t vb = COLOR_RGB565_TO_B5(img_pixel);

                            uint32_t malpha = 256 - alpha;
                            vr = ((vr * malpha) + (fr * alpha)) >> 8;
                            vg = ((vg * malpha) + (fg * alpha)) >> 8;
                            vb = ((vb * malpha) + (fb * alpha)) >> 8;

                            IMAGE_PUT_RGB565_PIXEL_FAST((uint16_t*)row_ptr, x_pixel, COLOR_R5_G6_B5_TO_RGB565(vr, vg, vb));
                        }
                        break;
                    }
                }
            }
            /*else if (intersect) {
                // skip empty pixels
                int delta = floorf(intersect->x) - x - 1;
                x += delta;
            }*/
        }
    }

    return xMax - xMin + 1;
}

void image_draw_ttf(image_t *img, const uint8_t* font, const char* text, uint32_t color, int x_off, int y_off, float size_pixels, int align, int valign, int justify) {
    ttf_font_instance ttf;
    ttf_init(&ttf, font);

    int len = strlen(text);
    int x_units = 0, y_units = 0;
    uint16_t units_per_em = __REV16(ttf.head->unitsPerEm);
    // float pixels_per_unit = size_pixels / units_per_em;

    if (valign >= 0) {
        uint32_t line_count = 1;

        // calculate how many lines in text
        for(int i=0; i < len; i++) if (text[i] == '\n') line_count++;

        if (valign == 1) {
            y_off = img->h - 1 - y_off;
            y_units -= line_count * units_per_em;
        }
        if (valign == 0) {
            y_off = (img->h >> 1) - y_off;
            y_units -= (line_count * units_per_em) >> 1;
        }
    }

    int32_t max_line_width = 0, line_width = 0;

    // Calculate maximum line width
    for(int i = 0; i < len; i++) {
        if (text[i] == '\n') {
            max_line_width = IM_MAX(max_line_width, line_width);
            line_width = 0;
        }
        else {
            uint16_t wch = (uint16_t)text[i] & 255;
            int glyph_index = ttf_table_cmap_char_to_index(&ttf, wch);
            longHorMetric metrics = ttf_get_horizontal_metrics(&ttf, glyph_index);
            line_width += __REV16(metrics.advanceWidth);
        }
    }
    max_line_width = IM_MAX(max_line_width, line_width);

    int x_off_i = x_off;

    for(int i=0; i < len; i++) {
        char ch = text[i];
        // Go to next line or is first line?
        if (ch == '\n' || i == 0 ) {

            // Calculate next lines width
            int32_t line_width = 0;
            for(int j = i + (ch == '\n' ? 1 : 0); j < len && text[j] != '\n'; j++) {
                uint16_t wch = (uint16_t)text[j] & 255;
                int glyph_index = ttf_table_cmap_char_to_index(&ttf, wch);
                longHorMetric metrics = ttf_get_horizontal_metrics(&ttf, glyph_index);
                line_width += __REV16(metrics.advanceWidth);
            }

            if (align == 1) {
                x_off = img->w - 1 - x_off_i;
                if (justify == 1) x_units = -line_width;
                else if (justify == 0) x_units = -((max_line_width + line_width) >> 1);
                else x_units = -max_line_width;
            } else if (align == 0) {
                x_off = (img->w >> 1) - x_off_i;
                if (justify == 1) x_units = (max_line_width >> 1) - line_width;
                else if (justify == 0) x_units = -(line_width >> 1);
                else x_units = -(max_line_width >> 1);
            }
            else {
                x_off = x_off_i;
                if (justify == 1) x_units = (max_line_width - line_width);
                else if (justify == 0) x_units = ((max_line_width - line_width) >> 1);
                else x_units = 0;
            }

            if (ch == '\n') {
                y_units += units_per_em;
            }
        }

        if (ch != '\n') {
            // Draw glyph
            int glyph_index = ttf_table_cmap_char_to_index(&ttf, (uint16_t)ch & 255);

            longHorMetric metrics = ttf_get_horizontal_metrics(&ttf, glyph_index);
            
            fb_alloc_mark();
            ttf_draw_glyph(&ttf, img, glyph_index, color, x_off, y_off, size_pixels, x_units, y_units);
            fb_alloc_free_till_mark();

            int32_t advance_width_units = __REV16(metrics.advanceWidth);
            x_units += advance_width_units;
        }
    }
}

// Validate font on load
// TODO Pre-convert font outline
// TODO sort lines by min y-start, min y-end
// 
// TODO load fonts from an sd card and cache
// TODO render horizontally and vertically at low res/oversample?
// TODO grid fit
// TODO image hints relative center top left right
// TODO unlimited glyph coordinates
// TODO glyph shadow for readability
// TODO glyph border
// Ensure at least one pixel gap between glyphs
/*int32_t advance_width_pixels = advance_width_units * pixels_per_unit;
if (advance_width_pixels <= glyph_width_pixels) {
advance_width_units = (glyph_width_pixels + 1) / pixels_per_unit;
}*/