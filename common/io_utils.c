//
// Created by wentao on 2024/2/4.
//

#include "io_utils.h"

size_t file_size(FILE* fp) {
    fseek(fp, 0, SEEK_END);
    size_t size_bytes = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    return size_bytes;
}