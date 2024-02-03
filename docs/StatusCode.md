# Status Code
- A 7-digit code:  showing the status of a variant.
- Reflecting the reliability of the statistics.
- Design principals:
    - Traceable
    - Higher value -> higher uncertainty 

| Digit | Description               |
| ----- | ------------------------- |
| 1,2   | Genome_build              |
| 3     | rsID & SNPID              |
| 4     | CHR, POS                  |
| 5     | EA, NEA                   |
| 6     | REF-NEA Alignment         |
| 7     | Palindromic SNPs + Indels |

## Summary

```
mysumstats.summary()
```

![image](https://user-images.githubusercontent.com/40289485/211861728-bc56389d-84a9-4929-95ab-0517ac063dfd.png)


## Look up the status code

```
mysumstats.lookup_status()
```

![image](https://user-images.githubusercontent.com/40289485/211861846-8309be9f-05ea-456e-ad8a-3265872826f9.png)


## Reference table
![image](https://user-images.githubusercontent.com/40289485/196681586-eb79707c-d866-4393-a0f5-b825686c9b04.png)


