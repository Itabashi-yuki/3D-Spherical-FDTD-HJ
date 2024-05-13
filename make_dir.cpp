#include <iostream>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <sstream>
#include "fdtd3d.h"

// namespace fs = std::filesystem;

std::string global_dirName;

void make_dir(){
    // 現在時刻を取得
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);

    // フォーマットを指定して時刻を文字列に変換
    std::tm bt = *std::localtime(&time_t_now);
    std::ostringstream oss;
    oss << std::put_time(&bt, "%Y%m%d_%H%M");

    // ディレクトリ名として使用
    std::string dirName;
    dirName = oss.str();

    // ディレクトリ作成
    std::filesystem::create_directory("./data/" + dirName);

    global_dirName = dirName;

    // std::cout << "Created directory: " << dirName << std::endl;
}
