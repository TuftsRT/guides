---
tags: nlp
---

# Deep Learning

Most modern NLP are language models. While a language model can be anything that predicts new text, in practice, these models are based on deep neural networks. In these notebooks, we will start from the very basics of training a language model from scratch and build up to training much larger and more advanced architectures. Importantly, this material is very hard and not intended, as of now, to be comprehesive. Instead, I recommend you work through these notebooks alongside other deep learning resources. For example, Jeremy Howard's [Practical Deep Learning for Coders](https://www.youtube.com/watch?v=8SF_h3xF3cE&list=PLfYUBJiXbdtSvpQjSnJJ_PmDQB_VyT5iU&ab_channel=JeremyHoward) is incredibly helpful. The bottom line is that no one was born knowing deep learning and that it can help to see this material presented in several different ways.

```{gallery-grid}
---
grid-columns: 1
---
- header: "{fas}`book` Training Named Entity Recognition models with `spaCy` "
  content: "In previous notebooks, we've seen how to use `spaCy` for NER, but how did those pretrained models get trained in the first place? This notebook shows you how!"
  link: "ner-spacy.html"

- header: "{fas}`book` An Introduction to Word Vectors with `gensim`"
  content: "Word vectors are the foundation of neual language models. The word2vec algorithm (Mikolov et al., 2013) pioneered this approach and the `gensim` package helps us implement it."
  link: "w2v-gensim.html"

- header: "{fas}`book` word2vec from Scratch "
  content: "Learn how to implement the word2vec algorithm in `numpy`. Start your foundation in training your own neural nets with this very simple architecture."
  link: "w2v-from-scratch.html"

- header: "{fas}`book` Building a machine translator with `pytorch` "
  content: "Get started with using the popluar neural network frameworking `pytorch`. Thisof the other side of the transformer architecture, the encoder, pivotal in BERT-style models for information retrieval and organization. notebook was adapted from [this source](https://pytorch.org/tutorials/intermediate/seq2seq_translation_tutorial.html)."
  link: "rnn-pytorch.html"

- header: "{fas}`book` Decoder-only Transformer Models from Scratch"
  content: "Understanding the decoder-only transformer, the architecture behind GPTs and LLMs, is vital to navigating the contemporary machine learning and artificial intelligence. Work through this notebook to learn the basics. This notebook was adapted from this [this source](https://www.youtube.com/watch?v=kCc8FmEb1nY&ab_channel=AndrejKarpathy)."
  link: "decoder-pytorch.html"

- header: "{fas}`book` Encoder-only Transformer Models for Scratch"
  content: "Develop your knowledge of the other side of the transformer architecture, the encoder, pivotal in BERT-style models for information retrieval and organization."
  link: "encoder-pytorch.html"

- header: "{fas}`book` The Byte-Pair Tokenization Algorithm"
  content: "Alongside model develop is text tokenization. In this notebook, learn modern techniques on how to build a custom tokenizer."
  link: "tokenizers.html"

```
