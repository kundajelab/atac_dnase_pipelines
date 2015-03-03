private int get_line_cnt(File f) { // count number of lines in a text file
    int i = 0
    f.eachLine { i++ }
    return i
}
