
#define MIN(a,b) (((a)<(b))?(a):(b))
typedef struct dictionary item;
static __inline item* push(unsigned int key, item* curr){
  item* head;
  head = (item*)malloc(sizeof(item));
  head->key = key;
  head->value = 0;
  head->next = curr;
  return head;
}
static __inline item* find(item* head, unsigned int key){
  item* iterator = head;
  while(iterator){
    if(iterator->key == key){
      return iterator;
    }
    iterator = iterator->next;
  }
  return NULL;
}
static __inline item* uniquePush(item* head, unsigned int key){
  item* iterator = head;
  while(iterator){
    if(iterator->key == key){
      return head;
    }
    iterator = iterator->next;
  }
  return push(key, head);
}
static void dict_free(item* head){
  item* iterator = head;
  while(iterator){
    item* temp = iterator;
    iterator = iterator->next;
    free(temp);
  }
  head = NULL;
}

static int dl_dist(int64_t src_input, int64_t tgt_input, unsigned int maxDistance){
  
  const unsigned int x = 16; // Hardcoded sequence length
  const unsigned int y = 16; // Hardcoded sequence length
  
  std::vector<int64_t> src_input_vec = {src_input}; // Wrap src_input in a vector
  std::vector<int64_t> tgt_input_vec = {tgt_input}; // Wrap tgt_input in a vector
  
  std::vector<unsigned int> src_vec = bits_to_uint_cpp(src_input_vec, x);
  std::vector<unsigned int> tgt_vec = bits_to_uint_cpp(tgt_input_vec, y);
  
  unsigned int* src = src_vec.data();  // Equivalent to &src_vec[0]
  unsigned int* tgt = tgt_vec.data();  // Equivalent to &tgt_vec[0]
  
  item *head = NULL;
  unsigned int swapCount, swapScore, targetCharCount, i, j;
  unsigned int *scores = (unsigned int*)malloc( (x + 2) * (y + 2) * sizeof(unsigned int) );
  unsigned int score_ceil = x + y;
  unsigned int curr_score;
  unsigned int diff = x > y ? (x - y) : (y - x);
  if(maxDistance != 0 && diff > maxDistance) {
    free(scores);
    return -1;
  }
  scores[0] = score_ceil;
  scores[1 * (y + 2) + 0] = score_ceil;
  scores[0 * (y + 2) + 1] = score_ceil;
  scores[1 * (y + 2) + 1] = 0;
  head = uniquePush(uniquePush(head, src[0]), tgt[0]);
  for(i = 1; i <= x; i++) {
    if(i < x)
      head = uniquePush(head, src[i]);
    scores[(i + 1) * (y + 2) + 1] = i;
    scores[(i + 1) * (y + 2) + 0] = score_ceil;
    swapCount = 0;
    for(j = 1; j <= y; j++) {
      if(i == 1) {
        if(j < y)
          head = uniquePush(head, tgt[j]);
        scores[1 * (y + 2) + (j + 1)] = j;
        scores[0 * (y + 2) + (j + 1)] = score_ceil;
      }
      curr_score = 0;
      targetCharCount = find(head, tgt[j - 1])->value;
      swapScore = scores[targetCharCount * (y + 2) + swapCount] + i - targetCharCount - 1 + j - swapCount;
      if(src[i - 1] != tgt[j - 1]){
        scores[(i + 1) * (y + 2) + (j + 1)] = MIN(swapScore, (MIN(scores[i * (y + 2) + j], MIN(scores[(i + 1) * (y + 2) + j], scores[i * (y + 2) + (j + 1)])) + 1));
      } else {
        swapCount = j;
        scores[(i + 1) * (y + 2) + (j + 1)] = MIN(scores[i * (y + 2) + j], swapScore);
      }
      curr_score = MIN(curr_score, scores[(i + 1) * (y + 2) + (j + 1)]);
    }
    if(maxDistance != 0 && curr_score > maxDistance) {
      dict_free(head);
      free(scores);
      return -1;
    }
    find(head, src[i - 1])->value = i;
  }
  {
    unsigned int score = scores[(x + 1) * (y + 2) + (y + 1)];
    dict_free(head);
    free(scores);
    return (maxDistance != 0 && maxDistance < score) ? (-1) : score;
  }
}