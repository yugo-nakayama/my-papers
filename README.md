# Papers Repository (Monorepo)

複数本の論文・プレプリント（arXiv投稿原稿等）を一元管理するモノレポ構成リポジトリです。各論文ディレクトリが自己完結し、arXiv向けのフラット構造パッケージ作成に対応しています。

## ディレクトリ構成

```text
my-papers/
├── .gitignore
├── LICENSE
├── README.md                # 本ファイル（論文一覧と概要）
├── scripts/                 # 共通ビルド・ユーティリティスクリプト
│   └── package_arxiv.sh    # arXiv提出用フラット化・パッケージングスクリプト
└── paper-<year>-<slug>/     # 個別論文ディレクトリ
    ├── main.tex             # メインTeXファイル
    ├── references.bib       # 参考文献（BibTeX）
    ├── figures/             # 図版ファイル
    ├── style/               # cls, styファイル等
    ├── Makefile             # ビルド・パッケージ作成ターゲット
    └── README.md            # 個別論文のメモ・進捗・リンク
```

## 収録論文一覧

| ディレクトリ | タイトル / トピック | ステータス | arXiv / 公開先 |
|---|---|---|---|
| [`paper-sse-svm`](./paper-sse-svm) | Asymptotic Properties of SVM with SSE Model (SSE_SVM) | 執筆 / 公開中 | [yugo-nakayama/SSE_SVM](https://github.com/yugo-nakayama/SSE_SVM) |

## 使い方

### 1. 新しい論文の作成
`paper-template` ディレクトリをコピーして利用してください:
```bash
cp -r paper-template paper-YYYY-topic-name
```

### 2. コンパイル
各論文ディレクトリに移動し、`make` を実行します:
```bash
cd paper-YYYY-topic-name
make
```

### 3. arXiv提出用アーカイブの作成
```bash
make arxiv
```
`arxiv_submission.tar.gz` が生成されます。arXivの投稿フォームにそのままアップロード可能です。

