<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Protocols</title>
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            line-height: 1.6;
            padding: 2rem;
            max-width: 800px;
            margin: 0 auto;
        }
        h1 {
            border-bottom: 2px solid #eee;
            padding-bottom: 0.5rem;
            margin-top: 2rem;
        }
        a {
            text-decoration: none;
            color: #007bff;
            transition: color 0.3s;
        }
        a:hover {
            color: #0056b3;
        }
    </style>
</head>
<body>

<header>
    <h1>Protocols</h1>
</header>

<main>
    <h1><a onclick="checkPassword()">Phylogenomics</a></h1>  
    <h1><a href="../protocols/transcriptome.html" target="_blank">Transcriptome</a></h1>  
    <h1><a href="../protocols/Molecular Convergence (conv_cal).html" target="_blank">Molecular Convergence (conv_cal)</a></h1>  
    <h1><a href="../protocols/Relaxed Selection Test.html" target="_blank">Relaxed Selection Test</a></h1>  
</main>

<script src="https://cdn.jsdelivr.net/npm/js-sha256@0.9.0/src/sha256.min.js"></script>

<script>
    function checkPassword() {
        const obfuscatedHashParts = [
            '1003f5', '9c15', '179cf021', '9817', 
            '6e9c', 'dd0333', '936c53e8', '3921ecde',
            '31003f56', 'e9c9c', '463e77c8', '3197ef'
        ];
        const getCorrectHash = () => {
            return obfuscatedHashParts[6] + obfuscatedHashParts[7] + obfuscatedHashParts[5] + 
                   '2463e77c8' + obfuscatedHashParts[8] + obfuscatedHashParts[4] + 
                   obfuscatedHashParts[9] + obfuscatedHashParts[2] + obfuscatedHashParts[11];
        };
        const correctHash = getCorrectHash();
        const userPassword = prompt("请输入访问密码：");
        if (userPassword === null || userPassword === "") {
            return;
        }
        const userHash = sha256(userPassword);
        if (userHash === correctHash) {
            window.location.href = "../protocols/Phylogenomics.html";
        } else {
            alert("密码错误，请重试！");
        }
    }
</script>

</body>
</html>
